import os
import subprocess
import threading
import wx

import preferences as pref

class MapReceptorFrame(wx.Frame):
    def __init__(self):
        super().__init__(parent=None, title="Map Receptor (AutoGrid)", size=(650, 480))
        self.worker = None
        self.init_ui()
        self.Bind(wx.EVT_CLOSE, self.on_close)
        self.Centre()

    def on_close(self, event):
        if self.worker is not None and self.worker.is_alive():
            wx.MessageBox(
                "Wait for the active mapping to finish before closing this window.",
                "Mapping in progress",
                wx.OK | wx.ICON_INFORMATION,
            )
            event.Veto()
            return
        event.Skip()

    def init_ui(self):
        panel = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        # --- Receptor GPF File Selection ---
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        file_label = wx.StaticText(panel, label="Receptor GPF File:")
        self.file_picker = wx.FilePickerCtrl(
            panel,
            path=os.getcwd(),
            message="Select receptor GPF file",
            wildcard="AutoGrid parameter files (*.gpf)|*.gpf|All files (*.*)|*.*",
        )
        hbox1.Add(file_label, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=8)
        hbox1.Add(self.file_picker, proportion=1, flag=wx.EXPAND)
        vbox.Add(hbox1, flag=wx.EXPAND | wx.ALL, border=10)

        # --- Run Button ---
        self.run_btn = wx.Button(panel, label="Run AutoGrid")
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
        gpf_file = self.file_picker.GetPath().strip()

        if not gpf_file or not os.path.isfile(gpf_file):
            wx.MessageBox("Please select an existing receptor GPF file.", "Error", wx.OK | wx.ICON_ERROR)
            return

        autogrid = os.path.join(str(pref.AutoDock_dir), 'AutoGrid', 'Linux64', 'autogrid4')
        if not os.path.isfile(autogrid):
            wx.MessageBox(
                f"autogrid4 not found:\n{autogrid}\n\nCheck AutoDock_dir in preferences.py.",
                "Error", wx.OK | wx.ICON_ERROR)
            return

        # Disable controls during execution
        self.run_btn.Disable()
        self.gauge.SetValue(0)
        self.log_txt.Clear()

        # Run mapping in a separate thread to keep the GUI active
        thread = threading.Thread(target=self.run_mapping, args=(gpf_file, autogrid))
        # Keep the process alive until AutoGrid has finished writing maps.
        thread.daemon = False
        self.worker = thread
        thread.start()

    def run_mapping(self, gpf_file, autogrid):
        dir_name = os.path.dirname(os.path.abspath(gpf_file))
        base_name = os.path.splitext(os.path.basename(gpf_file))[0]
        glg_file = os.path.join(dir_name, f"{base_name}.glg")

        cmd = [
            autogrid,
            '-p', gpf_file,
            '-l', glg_file,
        ]

        self.log(f"Mapping receptor: {gpf_file}")
        self.log("Running: " + " ".join(cmd) + "\n" + "-"*50)

        try:
            process = subprocess.Popen(
                cmd,
                cwd=dir_name,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True
            )
            stdout, _ = process.communicate()

            if stdout.strip():
                self.log(stdout.strip())

            if process.returncode == 0 and os.path.isfile(glg_file):
                wx.CallAfter(self.gauge.SetValue, 100)
                self.log(f"    --> Successfully created: {glg_file}")
                self.log("\n" + "="*50 + "\nReceptor mapping finished.")
            else:
                self.log(f"\n    --> ERROR: autogrid4 exited with code {process.returncode}")

        except Exception as e:
            self.log(f"    --> Failed to execute autogrid4: {str(e)}")

        wx.CallAfter(self.run_btn.Enable)


if __name__ == "__main__":
    app = wx.App(False)
    frame = MapReceptorFrame()
    frame.Show()
    app.MainLoop()
