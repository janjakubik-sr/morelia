import os
import glob
import re
import threading
import wx

class Smi2PdbqtFrame(wx.Frame):
    def __init__(self):
        super().__init__(parent=None, title="SMILES to PDBQT Converter (RDKit + Meeko)", size=(650, 480))
        self.worker = None
        self.init_ui()
        self.Bind(wx.EVT_CLOSE, self.on_close)
        self.Centre()

    def on_close(self, event):
        if self.worker is not None and self.worker.is_alive():
            wx.MessageBox(
                "Wait for the active conversion to finish before closing this window.",
                "Conversion in progress",
                wx.OK | wx.ICON_INFORMATION,
            )
            event.Veto()
            return
        event.Skip()

    def init_ui(self):
        panel = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        # --- Directory Selection ---
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        dir_label = wx.StaticText(panel, label="Working Directory:")
        self.dir_picker = wx.DirPickerCtrl(panel, path=os.getcwd(), message="Select search folder")
        hbox1.Add(dir_label, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=8)
        hbox1.Add(self.dir_picker, proportion=1, flag=wx.EXPAND)
        vbox.Add(hbox1, flag=wx.EXPAND | wx.ALL, border=10)

        # --- Glob Pattern Input ---
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        pattern_label = wx.StaticText(panel, label="File Glob Pattern:")
        self.pattern_txt = wx.TextCtrl(panel, value="*.smi")
        hbox2.Add(pattern_label, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=8)
        hbox2.Add(self.pattern_txt, proportion=1, flag=wx.EXPAND)
        vbox.Add(hbox2, flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, border=10)

        # --- Convert Button ---
        self.run_btn = wx.Button(panel, label="Start Conversion")
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
        work_dir = self.dir_picker.GetPath()
        pattern = self.pattern_txt.GetValue().strip()

        if not pattern:
            wx.MessageBox("Please enter a valid file glob pattern.", "Error", wx.OK | wx.ICON_ERROR)
            return

        # Disable controls during execution
        self.run_btn.Disable()
        self.gauge.SetValue(0)
        self.log_txt.Clear()

        # Run conversion in a separate thread to keep the GUI active
        thread = threading.Thread(target=self.run_conversion, args=(work_dir, pattern))
        # Keep the process alive until all ligands have been written.
        thread.daemon = False
        self.worker = thread
        thread.start()

    def run_conversion(self, work_dir, pattern):
        # Import chemistry backends lazily so the GUI can start and report
        # missing dependencies instead of crashing on launch.
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            from meeko import MoleculePreparation
        except ImportError as e:
            self.log(f"ERROR: Required Python module not found: {e.name}")
            self.log("Please install 'rdkit' and 'meeko' in the active Python environment.")
            wx.CallAfter(self.run_btn.Enable)
            return

        search_path = os.path.join(work_dir, pattern)
        files = glob.glob(search_path)

        if not files:
            self.log(f"No files matched pattern: {pattern}")
            wx.CallAfter(self.run_btn.Enable)
            return

        total_files = len(files)
        self.log(f"Found {total_files} file(s) matching pattern.\n" + "-"*50)

        for idx, file_path in enumerate(files, start=1):
            dir_name = os.path.dirname(file_path)
            base_name = os.path.splitext(os.path.basename(file_path))[0]

            self.log(f"[{idx}/{total_files}] Processing: {file_path}")

            try:
                with open(file_path, "r") as smi_file:
                    lines = smi_file.readlines()
            except Exception as e:
                self.log(f"    --> ERROR reading {file_path}: {str(e)}")
                wx.CallAfter(self.gauge.SetValue, int((idx / total_files) * 100))
                continue

            mol_count = 0
            for line_num, line in enumerate(lines, start=1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                # Standard .smi format: SMILES optionally followed by a name
                fields = line.split()
                smiles = fields[0]
                mol_name = fields[1] if len(fields) > 1 else f"{base_name}_{line_num}"
                # Make the name safe for use as a file name
                mol_name = re.sub(r'[^A-Za-z0-9_.-]', '_', mol_name)
                out_pdbqt = os.path.join(dir_name, f"{mol_name}.pdbqt")

                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    self.log(f"    --> SKIP line {line_num}: invalid SMILES '{smiles}'")
                    continue

                try:
                    mol = Chem.AddHs(mol)
                    if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0:
                        self.log(f"    --> SKIP line {line_num}: 3D embedding failed for '{smiles}'")
                        continue
                    AllChem.MMFFOptimizeMolecule(mol)

                    preparer = MoleculePreparation()
                    preparer.prepare(mol)
                    pdbqt_string = preparer.write_pdbqt_string()

                    tmp_pdbqt = out_pdbqt + ".part"
                    with open(tmp_pdbqt, "w") as f:
                        f.write(pdbqt_string)
                    os.replace(tmp_pdbqt, out_pdbqt)

                    mol_count += 1
                    self.log(f"    --> Created: {out_pdbqt}")
                except Exception as e:
                    self.log(f"    --> ERROR on line {line_num} ('{smiles}'): {str(e)}")

            if mol_count == 0:
                self.log(f"    --> WARNING: No ligands converted from {file_path}")

            # Update Progress Bar
            progress = int((idx / total_files) * 100)
            wx.CallAfter(self.gauge.SetValue, progress)

        self.log("\n" + "="*50 + "\nAll conversion tasks finished.")
        wx.CallAfter(self.run_btn.Enable)


if __name__ == "__main__":
    app = wx.App(False)
    frame = Smi2PdbqtFrame()
    frame.Show()
    app.MainLoop()
