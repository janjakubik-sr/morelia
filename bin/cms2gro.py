import os
import glob
import subprocess
import threading
import wx

class VMDConverterFrame(wx.Frame):
    def __init__(self):
        super().__init__(parent=None, title="Desmond .cms to GROMACS .gro Converter", size=(650, 480))
        self.init_ui()
        self.Centre()

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
        thread.daemon = True
        thread.start()

    def run_conversion(self, work_dir, pattern):
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
            out_gro = os.path.join(dir_name, f"{base_name}.gro")

            self.log(f"[{idx}/{total_files}] Converting: {file_path}")

            # Construct Tcl script sent directly into VMD's stdin
            vmd_tcl = f"""mol new "{file_path}" type mae waitfor all
set sel [atomselect top "all"]
$sel writegro "{out_gro}"
quit
"""

            try:
                process = subprocess.Popen(
                    ["vmd", "-dispdev", "text"],
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True
                )
                stdout, stderr = process.communicate(input=vmd_tcl)

                if process.returncode == 0:
                    self.log(f"    --> Successfully created: {out_gro}")
                else:
                    self.log(f"    --> ERROR processing {file_path}:\n{stderr}")

            except Exception as e:
                self.log(f"    --> Failed to execute VMD: {str(e)}")

            # Update Progress Bar
            progress = int((idx / total_files) * 100)
            wx.CallAfter(self.gauge.SetValue, progress)

        self.log("\n" + "="*50 + "\nAll conversion tasks finished.")
        wx.CallAfter(self.run_btn.Enable)


if __name__ == "__main__":
    app = wx.App(False)
    frame = VMDConverterFrame()
    frame.Show()
    app.MainLoop()