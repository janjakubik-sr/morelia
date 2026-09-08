import os
import csv
import glob
import re
import subprocess
import threading
import wx

import preferences as pref

FLOAT_RE = r'[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?'

# AD4 atom type -> element symbol for PDB element columns
AD4_METALS = {'FE': 'Fe', 'MG': 'Mg', 'ZN': 'Zn', 'CA': 'Ca', 'MN': 'Mn', 'CU': 'Cu'}

def adtype_to_element(atom_type):
    t = atom_type.strip()
    if not t:
        return 'C'
    if t == 'A':            # aromatic carbon
        return 'C'
    if t in AD4_METALS:
        return AD4_METALS[t]
    if len(t) == 2 and t[0].isupper() and t[1].islower():
        return t            # Cl, Br, ...
    return t[0]             # C, N, NA, OA, SA, HD, F, H, P, I, S, ...

def pdbqt_to_pdb_lines(pdbqt_path, record=None, serial_offset=0):
    """Convert ATOM/HETATM records of a PDBQT file to PDB lines.

    Strips partial charge and AD4 atom type columns, writes the element
    symbol into columns 77-78, and optionally renumbers serials.
    Returns (lines, max_serial).
    """
    lines = []
    max_serial = 0
    with open(pdbqt_path, 'r') as f:
        for raw in f:
            if not (raw.startswith('ATOM') or raw.startswith('HETATM')):
                continue
            tokens = raw.split()
            element = adtype_to_element(tokens[-1]) if tokens else 'C'
            base = raw.rstrip('\n')[:66].ljust(76)
            if record is not None:
                base = record.ljust(6) + base[6:]
            if serial_offset:
                serial = int(base[6:11]) + serial_offset
                base = base[:6] + str(serial % 100000).rjust(5) + base[11:]
            max_serial = max(max_serial, int(base[6:11]))
            lines.append(base + element.rjust(2))
    return lines, max_serial

def parse_dlg(dlg_path):
    """Parse an AutoDock DLG file.

    Returns (clusters, poses) where clusters are dicts with keys
    rank, lowest_energy, run, mean_energy, num_in_cluster and poses are
    dicts with keys rank, sub_rank, run, energy, cluster_rmsd, ref_rmsd.
    """
    cluster_re = re.compile(
        r'^\s*(\d+)\s*\|\s*(%s)\s*\|\s*(\d+)\s*\|\s*(%s)\s*\|\s*(\d+)\s*\|'
        % (FLOAT_RE, FLOAT_RE))
    pose_re = re.compile(
        r'^\s*(\d+)\s+(\d+)\s+(\d+)\s+(%s)\s+(%s)\s+(%s)\s+RANKING\s*$'
        % (FLOAT_RE, FLOAT_RE, FLOAT_RE))
    clusters, poses = [], []
    with open(dlg_path, 'r', errors='replace') as f:
        for line in f:
            m = cluster_re.match(line)
            if m:
                clusters.append({
                    'rank': int(m.group(1)),
                    'lowest_energy': float(m.group(2)),
                    'run': int(m.group(3)),
                    'mean_energy': float(m.group(4)),
                    'num_in_cluster': int(m.group(5)),
                })
                continue
            m = pose_re.match(line)
            if m:
                poses.append({
                    'rank': int(m.group(1)),
                    'sub_rank': int(m.group(2)),
                    'run': int(m.group(3)),
                    'energy': float(m.group(4)),
                    'cluster_rmsd': float(m.group(5)),
                    'ref_rmsd': float(m.group(6)),
                })
    return clusters, poses

def receptor_pdbqt_from_fld(fld_path):
    """Find the receptor PDBQT referenced in the .fld header."""
    fld_dir = os.path.dirname(fld_path)
    try:
        with open(fld_path, 'r', errors='replace') as f:
            for line in f:
                if line.startswith('#MACROMOLECULE'):
                    candidate = line.split(None, 1)[1].strip()
                    path = candidate if os.path.isabs(candidate) \
                        else os.path.join(fld_dir, candidate)
                    if os.path.isfile(path):
                        return path
    except Exception:
        pass
    # Conventional fallback: <base>.maps.fld -> <base>.pdbqt
    base = os.path.basename(fld_path)
    for suffix in ('.maps.fld', '.fld'):
        if base.endswith(suffix):
            path = os.path.join(fld_dir, base[:-len(suffix)] + '.pdbqt')
            if os.path.isfile(path):
                return path
    return None

class BatchDockFrame(wx.Frame):
    def __init__(self):
        super().__init__(parent=None, title="Batch Docking (AutoDock-GPU)", size=(650, 520))
        self.worker = None
        self.init_ui()
        self.Bind(wx.EVT_CLOSE, self.on_close)
        self.Centre()

    def on_close(self, event):
        if self.worker is not None and self.worker.is_alive():
            wx.MessageBox(
                "Wait for the active docking to finish before closing this window.",
                "Docking in progress",
                wx.OK | wx.ICON_INFORMATION,
            )
            event.Veto()
            return
        event.Skip()

    def init_ui(self):
        panel = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        # --- Receptor FLD File Selection ---
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        file_label = wx.StaticText(panel, label="Receptor FLD File:")
        self.file_picker = wx.FilePickerCtrl(
            panel,
            path=os.getcwd(),
            message="Select receptor grid maps .fld file",
            wildcard="AutoDock map descriptor files (*.fld)|*.fld|All files (*.*)|*.*",
        )
        hbox1.Add(file_label, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=8)
        hbox1.Add(self.file_picker, proportion=1, flag=wx.EXPAND)
        vbox.Add(hbox1, flag=wx.EXPAND | wx.ALL, border=10)

        # --- Ligand Glob Pattern Input ---
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        pattern_label = wx.StaticText(panel, label="Ligand Glob Pattern:")
        self.pattern_txt = wx.TextCtrl(panel, value="*/*.pdbqt")
        hbox2.Add(pattern_label, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=8)
        hbox2.Add(self.pattern_txt, proportion=1, flag=wx.EXPAND)
        vbox.Add(hbox2, flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, border=10)

        # --- Number of Runs Input ---
        hbox3 = wx.BoxSizer(wx.HORIZONTAL)
        hbox3.Add(wx.StaticText(panel, label="Number of Runs:"),
                  flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=8)
        self.nrun_txt = wx.TextCtrl(panel, value="888", size=(70, -1))
        hbox3.Add(self.nrun_txt, flag=wx.ALIGN_CENTER_VERTICAL)
        vbox.Add(hbox3, flag=wx.LEFT | wx.RIGHT | wx.BOTTOM, border=10)

        # --- Run Button ---
        self.run_btn = wx.Button(panel, label="Start Docking")
        self.run_btn.Bind(wx.EVT_BUTTON, self.on_start)
        vbox.Add(self.run_btn, flag=wx.ALIGN_CENTER | wx.BOTTOM, border=10)

        # --- Progress & Logging ---
        self.gauge = wx.Gauge(panel, range=100)
        vbox.Add(self.gauge, flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, border=10)

        log_label = wx.StaticText(panel, label="Console Output:")
        vbox.Add(log_label, flag=wx.LEFT, border=10)

        self.log_txt = wx.TextCtrl(panel, style=wx.TE_MULTILINE | wx.TE_READONLY)
        vbox.Add(self.log_txt, proportion=1, flag=wx.EXPAND | wx.ALL, border=10)

        panel.SetSizer(vbox)

    def log(self, message):
        """Append log text smoothly in the main thread."""
        wx.CallAfter(self.log_txt.AppendText, message + "\n")

    def on_start(self, event):
        fld_file = self.file_picker.GetPath().strip()
        pattern = self.pattern_txt.GetValue().strip()

        if not fld_file or not os.path.isfile(fld_file):
            wx.MessageBox("Please select an existing receptor FLD file.", "Error", wx.OK | wx.ICON_ERROR)
            return
        if not pattern:
            wx.MessageBox("Please enter a valid ligand glob pattern.", "Error", wx.OK | wx.ICON_ERROR)
            return
        try:
            nrun = int(self.nrun_txt.GetValue())
            if nrun < 1:
                raise ValueError
        except ValueError:
            wx.MessageBox("Number of runs must be a positive integer.", "Error", wx.OK | wx.ICON_ERROR)
            return

        autodock = os.path.join(str(pref.AutoDock_dir), 'AutoDock-GPU', 'bin', 'autodock_gpu_64wi')
        if not os.path.isfile(autodock):
            wx.MessageBox(
                f"autodock_gpu_64wi not found:\n{autodock}\n\nCheck AutoDock_dir in preferences.py.",
                "Error", wx.OK | wx.ICON_ERROR)
            return

        # Disable controls during execution
        self.run_btn.Disable()
        self.gauge.SetValue(0)
        self.log_txt.Clear()

        # Run docking in a separate thread to keep the GUI active
        thread = threading.Thread(target=self.run_docking, args=(fld_file, pattern, nrun, autodock))
        # Keep the process alive until all dockings have finished.
        thread.daemon = False
        self.worker = thread
        thread.start()

    def run_docking(self, fld_file, pattern, nrun, autodock):
        work_dir = os.path.dirname(os.path.abspath(fld_file))
        fld_file = os.path.abspath(fld_file)
        ligands = sorted(glob.glob(os.path.join(work_dir, pattern)))

        if not ligands:
            self.log(f"No ligands matched pattern: {pattern}")
            wx.CallAfter(self.run_btn.Enable)
            return

        receptor_pdbqt = receptor_pdbqt_from_fld(fld_file)
        if receptor_pdbqt is None:
            self.log("WARNING: Receptor PDBQT not found; complex export will be skipped.")
        else:
            self.log(f"Receptor PDBQT for complex export: {receptor_pdbqt}")

        total = len(ligands)
        self.log(f"Found {total} ligand(s) matching pattern.\n" + "-"*50)

        for idx, ligand in enumerate(ligands, start=1):
            lig_dir = os.path.dirname(os.path.abspath(ligand))
            base = os.path.splitext(os.path.basename(ligand))[0]
            lig_id = os.path.splitext(os.path.relpath(ligand, work_dir))[0]
            resnam = os.path.join(lig_dir, f"{base}_dock")
            dlg_path = resnam + ".dlg"
            gbest_path = resnam + "-best.pdbqt"

            self.log(f"[{idx}/{total}] Docking: {lig_id}  ({nrun} runs)")

            cmd = [
                autodock,
                '--lfile', ligand,
                '--ffile', fld_file,
                '--nrun', str(nrun),
                '--resnam', resnam,
                '--gbest', '1',
            ]

            try:
                process = subprocess.Popen(
                    cmd,
                    cwd=lig_dir,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True
                )
                stdout, _ = process.communicate()

                if process.returncode != 0 or not os.path.isfile(dlg_path):
                    self.log(f"    --> ERROR docking {lig_id}:\n{stdout.strip()}")
                    wx.CallAfter(self.gauge.SetValue, int((idx / total) * 100))
                    continue

                clusters, poses = parse_dlg(dlg_path)
                if not clusters or not poses:
                    self.log(f"    --> WARNING: No clusters/poses parsed from {dlg_path}")
                for row in clusters:
                    row['ligand'] = lig_id
                for row in poses:
                    row['ligand'] = lig_id

                # Write sorted CSV tables into the ligand subfolder
                if clusters:
                    clusters.sort(key=lambda r: (r['lowest_energy'], r['rank']))
                    clusters_csv = os.path.join(lig_dir, f"{base}_clusters.csv")
                    with open(clusters_csv, 'w', newline='') as f:
                        writer = csv.writer(f)
                        writer.writerow(['ligand', 'cluster_rank', 'lowest_binding_energy',
                                         'mean_binding_energy', 'num_in_cluster', 'best_run'])
                        for r in clusters:
                            writer.writerow([r['ligand'], r['rank'], f"{r['lowest_energy']:.2f}",
                                             f"{r['mean_energy']:.2f}", r['num_in_cluster'], r['run']])
                    self.log(f"    --> Clusters table: {clusters_csv}")

                if poses:
                    poses.sort(key=lambda r: (r['energy'], r['rank'], r['sub_rank']))
                    poses_csv = os.path.join(lig_dir, f"{base}_poses.csv")
                    with open(poses_csv, 'w', newline='') as f:
                        writer = csv.writer(f)
                        writer.writerow(['ligand', 'cluster_rank', 'cluster_sub_rank', 'run',
                                         'binding_energy', 'cluster_rmsd', 'reference_rmsd'])
                        for r in poses:
                            writer.writerow([r['ligand'], r['rank'], r['sub_rank'], r['run'],
                                             f"{r['energy']:.2f}", f"{r['cluster_rmsd']:.2f}",
                                             f"{r['ref_rmsd']:.2f}"])
                    self.log(f"    --> Poses table: {poses_csv}")

                best = min(poses, key=lambda r: r['energy'])
                self.log(
                    f"    --> Best pose: {best['energy']:.2f} kcal/mol"
                    f" (run {best['run']}, cluster {best['rank']},"
                    f" cluster RMSD {best['cluster_rmsd']:.2f}, ref RMSD {best['ref_rmsd']:.2f})"
                )

                # Export top binding pose ligand-receptor complex to PDB
                if receptor_pdbqt is not None and os.path.isfile(gbest_path):
                    complex_pdb = os.path.join(lig_dir, f"{base}_complex.pdb")
                    rec_lines, rec_max = pdbqt_to_pdb_lines(receptor_pdbqt)
                    lig_lines, _ = pdbqt_to_pdb_lines(gbest_path, record='HETATM', serial_offset=rec_max)
                    tmp_pdb = complex_pdb + ".part"
                    with open(tmp_pdb, 'w') as f:
                        f.write("\n".join(rec_lines + lig_lines + ["END"]) + "\n")
                    os.replace(tmp_pdb, complex_pdb)
                    self.log(f"    --> Complex: {complex_pdb}")

            except Exception as e:
                self.log(f"    --> Failed to dock {lig_id}: {str(e)}")

            progress = int((idx / total) * 100)
            wx.CallAfter(self.gauge.SetValue, progress)

        self.log("\n" + "="*50 + "\nAll docking tasks finished.")
        wx.CallAfter(self.run_btn.Enable)


if __name__ == "__main__":
    app = wx.App(False)
    frame = BatchDockFrame()
    frame.Show()
    app.MainLoop()
