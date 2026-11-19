import os
import glob
import re
import shlex
import subprocess
import threading
import wx

# Conversions under this directory run on the file server "atlas" (sshfs
# mount) instead of locally, so the I/O-heavy VMD work stays server-side.
REMOTE_HOST = "atlas"
REMOTE_MOUNT = "/home/roshi/Data"


def is_remote_path(path):
    """True when *path* lives on the sshfs mount served by REMOTE_HOST."""
    abs_path = os.path.abspath(path)
    return abs_path == REMOTE_MOUNT or abs_path.startswith(REMOTE_MOUNT + os.sep)


def build_command(cmd_args, remote):
    """Return the command to run, wrapped in ssh when *remote* is true."""
    if remote:
        return ["ssh", REMOTE_HOST, " ".join(shlex.quote(arg) for arg in cmd_args)]
    return cmd_args


# .gro atom records are fixed format: 15 columns of residue/atom naming,
# a 5-column atom serial, then coordinates.
GRO_HEADER_COLS = 15
GRO_SERIAL_COLS = 5
GRO_COORD_COLS = 24  # three %8.3f fields; VMD's writegro emits no velocities


def wrap_gro_atom_serials(path):
    """Wrap overflowing .gro atom serials back into their 5-column field.

    VMD prints atom serials into the 5-column field without wrapping, so a
    system with more than 99999 atoms lets every serial from 100000 on
    spill into a sixth column.  That shifts the coordinate columns right
    and GROMACS aborts with "Something is wrong in the coordinate
    formatting".  GROMACS's own writer wraps serials modulo 100000
    (groio.cpp, write_hconf_indexed_p); the same convention is applied
    here.  Returns the number of atom records that were rewritten.
    """
    with open(path) as fh:
        lines = fh.readlines()

    natoms = int(lines[1])
    if len(lines) != natoms + 3:
        raise RuntimeError(
            f"{path}: expected {natoms} atom records, found {len(lines) - 3}"
        )

    fixed = 0
    for i in range(2, 2 + natoms):
        line = lines[i].rstrip("\n")
        tail = line[GRO_HEADER_COLS:]
        match = re.match(r"\s*(\d+)", tail)
        if match is None or len(tail) - match.end() != GRO_COORD_COLS:
            raise RuntimeError(
                f"{path}: unrecognized .gro atom record {i - 1}: {line!r}"
            )
        serial = int(match.group(1))
        if serial > 99999:
            lines[i] = (
                line[:GRO_HEADER_COLS]
                + f"{serial % 100000:{GRO_SERIAL_COLS}d}"
                + tail[match.end():]
                + "\n"
            )
            fixed += 1

    if fixed:
        with open(path, "w") as fh:
            fh.writelines(lines)
    return fixed


class VMDTrajectoryConverterFrame(wx.Frame):
    def __init__(self):
        super().__init__(parent=None, title="Desmond Trajectory to GROMACS .trr Converter", size=(670, 500))
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
        self.pattern_txt = wx.TextCtrl(panel, value="desmond_md_job_?/desmond_md_job_?.cms")
        hbox2.Add(pattern_label, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=8)
        hbox2.Add(self.pattern_txt, proportion=1, flag=wx.EXPAND)
        vbox.Add(hbox2, flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, border=10)

        # --- Frame Range Input ---
        frame_box = wx.BoxSizer(wx.HORIZONTAL)
        frame_box.Add(wx.StaticText(panel, label="Frame Range:"),
                      flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=8)

        self.beg_txt = wx.TextCtrl(panel, value="1", size=(70, -1))
        self.end_txt = wx.TextCtrl(panel, value="-1", size=(70, -1))
        self.skip_txt = wx.TextCtrl(panel, value="1", size=(70, -1))

        for label, control in (
            ("Beg:", self.beg_txt),
            ("End:", self.end_txt),
            ("Skip:", self.skip_txt),
        ):
            frame_box.Add(wx.StaticText(panel, label=label),
                          flag=wx.ALIGN_CENTER_VERTICAL | wx.LEFT | wx.RIGHT, border=5)
            frame_box.Add(control, flag=wx.ALIGN_CENTER_VERTICAL)

        vbox.Add(frame_box, flag=wx.LEFT | wx.RIGHT | wx.BOTTOM, border=10)

        # --- Convert Button ---
        self.run_btn = wx.Button(panel, label="Start Trajectory Conversion")
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

        try:
            beg = int(self.beg_txt.GetValue())
            end = int(self.end_txt.GetValue())
            skip = int(self.skip_txt.GetValue())
        except ValueError:
            wx.MessageBox("Beg, end, and skip must be integers.", "Error", wx.OK | wx.ICON_ERROR)
            return

        if beg < 0 or end < -1 or skip < 1:
            wx.MessageBox(
                "Beg must be 0 or greater, end must be -1 or greater, and skip must be 1 or greater.",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            return

        if end != -1 and end < beg:
            wx.MessageBox("End must be -1 or greater than or equal to beg.", "Error", wx.OK | wx.ICON_ERROR)
            return

        # Disable controls during execution
        self.run_btn.Disable()
        self.gauge.SetValue(0)
        self.log_txt.Clear()

        # Run conversion in a separate thread to keep GUI responsive
        thread = threading.Thread(
            target=self.run_conversion,
            args=(work_dir, pattern, beg, end, skip),
        )
        # Do not let closing the GUI terminate a conversion midway through a
        # trajectory file.
        thread.daemon = False
        self.worker = thread
        thread.start()

    def run_conversion(self, work_dir, pattern, beg, end, skip):
        search_path = os.path.join(work_dir, pattern)
        files = glob.glob(search_path)

        if not files:
            self.log(f"No files matched pattern: {pattern}")
            wx.CallAfter(self.run_btn.Enable)
            return

        total_files = len(files)
        self.log(f"Found {total_files} file(s) matching pattern.\n" + "-"*55)

        # Heavy trajectory I/O must happen where the files physically live.
        remote = is_remote_path(work_dir)
        if remote:
            self.log(f"Directory is on the sshfs mount; running VMD on {REMOTE_HOST}.")

        for idx, file_path in enumerate(files, start=1):
            dir_name = os.path.dirname(file_path)
            base_name = os.path.splitext(os.path.basename(file_path))[0]
            
            # Construct paths matching your bash script logic
            trj_dir = os.path.join(dir_name, f"{base_name}_trj")
            out_trr = os.path.join(dir_name, f"{base_name}.trr")
            out_gro = os.path.join(dir_name, f"{base_name}.gro")
            tmp_trr = out_trr + ".part"
            tmp_gro = out_gro + ".part"

            # Validate that trajectory folder exists before calling VMD
            if not os.path.exists(trj_dir):
                self.log(f"[{idx}/{total_files}] SKIP: Trajectory directory not found:\n    --> {trj_dir}")
                continue

            self.log(f"[{idx}/{total_files}] Converting trajectory: {base_name}_trj")

            # The CMS is VMD frame 0 and the DTR begins at frame 1.  Export the
            # GRO from the exact frame that becomes the first TRR frame.
            vmd_tcl = f"""mol new "{file_path}" type mae waitfor all
animate read dtr "{trj_dir}" waitfor all top
set sel [atomselect top "all"]
animate write trr "{tmp_trr}" beg {beg} end {end} skip {skip} sel $sel
animate goto {beg}
set gro_sel [atomselect top "all" frame {beg}]
$gro_sel writegro "{tmp_gro}"
quit
"""

            try:
                process = subprocess.Popen(
                    build_command(["vmd", "-dispdev", "text"], remote),
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True
                )
                stdout, stderr = process.communicate(input=vmd_tcl)

                if (
                    process.returncode == 0
                    and os.path.getsize(tmp_trr) > 0
                    and os.path.getsize(tmp_gro) > 0
                ):
                    fixed = wrap_gro_atom_serials(tmp_gro)
                    if fixed:
                        self.log(
                            f"    --> Wrapped {fixed} atom serials above the"
                            " 5-column .gro limit"
                        )
                    os.replace(tmp_trr, out_trr)
                    os.replace(tmp_gro, out_gro)
                    self.log(f"    --> Successfully created: {out_trr}")
                    self.log(f"    --> Synchronized frame {beg} structure: {out_gro}")
                else:
                    self.log(f"    --> ERROR processing {base_name}:\n{stderr}")

            except Exception as e:
                self.log(f"    --> Failed to execute VMD: {str(e)}")
            finally:
                for temporary_path in (tmp_trr, tmp_gro):
                    if os.path.exists(temporary_path):
                        os.remove(temporary_path)

            # Update Progress Bar
            progress = int((idx / total_files) * 100)
            wx.CallAfter(self.gauge.SetValue, progress)

        self.log("\n" + "="*55 + "\nAll trajectory conversion tasks finished.")
        wx.CallAfter(self.run_btn.Enable)


if __name__ == "__main__":
    app = wx.App(False)
    frame = VMDTrajectoryConverterFrame()
    frame.Show()
    app.MainLoop()
