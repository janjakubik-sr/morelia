import os
import glob
import re
import shlex
import subprocess
import threading
import wx

# Conversions under this directory run on the file server "atlas" (sshfs
# mount) instead of locally, so the I/O-heavy GROMACS work stays server-side.
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


def gmx_frame_count(path, remote=False):
    """Return the coordinate-frame count reported by ``gmx check``."""
    process = subprocess.run(
        build_command(["gmx", "check", "-f", path], remote),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    output = process.stdout + "\n" + process.stderr
    if process.returncode != 0:
        raise RuntimeError(f"gmx check failed for {path}:\n{output}")
    match = re.search(r"^\s*Coords\s+(\d+)\b", output, flags=re.MULTILINE)
    if match is None:
        raise RuntimeError(f"Could not determine frame count for {path} from gmx check")
    return int(match.group(1))

class GMXConvertConverterFrame(wx.Frame):
    def __init__(self):
        super().__init__(parent=None, title="GROMACS .trr to .xtc Compressor", size=(680, 530))
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
        self.pattern_txt = wx.TextCtrl(panel, value="desmond_md_job_?/desmond_md_job_?.trr")
        hbox2.Add(pattern_label, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=8)
        hbox2.Add(self.pattern_txt, proportion=1, flag=wx.EXPAND)
        vbox.Add(hbox2, flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, border=10)

        # --- Group Selection Input ---
        hbox3 = wx.BoxSizer(wx.HORIZONTAL)
        group_label = wx.StaticText(panel, label="GMX Group ID / Choice:")
        self.group_txt = wx.TextCtrl(panel, value="0")
        self.group_txt.SetToolTip("The group index or name GROMACS expects (e.g., 0 for System)")
        hbox3.Add(group_label, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=8)
        hbox3.Add(self.group_txt, proportion=1, flag=wx.EXPAND)
        vbox.Add(hbox3, flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, border=10)

        # --- Convert Button ---
        self.run_btn = wx.Button(panel, label="Start .xtc Compression")
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
        group_choice = self.group_txt.GetValue().strip()

        if not pattern:
            wx.MessageBox("Please enter a valid file glob pattern.", "Error", wx.OK | wx.ICON_ERROR)
            return

        # Disable controls during execution
        self.run_btn.Disable()
        self.gauge.SetValue(0)
        self.log_txt.Clear()

        # Run conversion in a separate thread to keep GUI responsive
        thread = threading.Thread(target=self.run_conversion, args=(work_dir, pattern, group_choice))
        # A daemon worker can be killed as the GUI exits, leaving a truncated
        # but superficially valid XTC behind.
        thread.daemon = False
        self.worker = thread
        thread.start()

    def run_conversion(self, work_dir, pattern, group_choice):
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
            self.log(f"Directory is on the sshfs mount; running GROMACS on {REMOTE_HOST}.")

        for idx, file_path in enumerate(files, start=1):
            dir_name = os.path.dirname(file_path)
            base_name = os.path.splitext(os.path.basename(file_path))[0]
            
            gro_file = os.path.join(dir_name, f"{base_name}.gro")
            out_xtc = os.path.join(dir_name, f"{base_name}.xtc")
            tmp_xtc = os.path.join(
                dir_name, f".{base_name}.{os.getpid()}.part.xtc"
            )

            # Check if matching structure (.gro) exists
            if not os.path.exists(gro_file):
                self.log(f"[{idx}/{total_files}] SKIP: Matching structure file not found:\n    --> {gro_file}")
                continue

            self.log(f"[{idx}/{total_files}] Compressing to .xtc: {base_name}.trr")

            # Command array for GROMACS trjconv
            cmd = [
                "gmx", "trjconv",
                "-f", file_path,
                "-s", gro_file,
                "-o", tmp_xtc
            ]

            try:
                # Pipe input to feed the group ID option automatically when GROMACS prompts
                process = subprocess.Popen(
                    build_command(cmd, remote),
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True
                )

                # Send the group choice followed by a newline and finish communication
                stdout, stderr = process.communicate(input=f"{group_choice}\n")

                if process.returncode == 0:
                    input_frames = gmx_frame_count(file_path, remote)
                    output_frames = gmx_frame_count(tmp_xtc, remote)
                    if output_frames != input_frames:
                        raise RuntimeError(
                            "Frame-count mismatch: "
                            f"TRR has {input_frames}, temporary XTC has {output_frames}"
                        )
                    os.replace(tmp_xtc, out_xtc)
                    self.log(
                        f"    --> Successfully created and validated "
                        f"{output_frames} frames: {out_xtc}"
                    )
                else:
                    self.log(f"    --> ERROR processing {base_name}:\n{stderr}")

            except Exception as e:
                self.log(f"    --> Failed to execute GROMACS: {str(e)}")
            finally:
                if os.path.exists(tmp_xtc):
                    os.remove(tmp_xtc)

            # Update Progress Bar
            progress = int((idx / total_files) * 100)
            wx.CallAfter(self.gauge.SetValue, progress)

        self.log("\n" + "="*55 + "\nAll .xtc compression tasks finished.")
        wx.CallAfter(self.run_btn.Enable)


if __name__ == "__main__":
    app = wx.App(False)
    frame = GMXConvertConverterFrame()
    frame.Show()
    app.MainLoop()
