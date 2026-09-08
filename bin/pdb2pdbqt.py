import os
import subprocess
import threading
import wx

import preferences as pref

class Pdb2PdbqtFrame(wx.Frame):
    def __init__(self):
        super().__init__(parent=None, title="Receptor PDB to PDBQT Converter (Meeko)", size=(650, 520))
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

        # --- Receptor PDB File Selection ---
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        file_label = wx.StaticText(panel, label="Receptor PDB File:")
        self.file_picker = wx.FilePickerCtrl(
            panel,
            path=os.getcwd(),
            message="Select receptor PDB file",
            wildcard="PDB files (*.pdb)|*.pdb|All files (*.*)|*.*",
        )
        hbox1.Add(file_label, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=8)
        hbox1.Add(self.file_picker, proportion=1, flag=wx.EXPAND)
        vbox.Add(hbox1, flag=wx.EXPAND | wx.ALL, border=10)

        # --- Box Center Input ---
        center_box = wx.BoxSizer(wx.HORIZONTAL)
        center_box.Add(wx.StaticText(panel, label="Box Center:"),
                       flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=8)
        self.cx_txt = wx.TextCtrl(panel, value="-2.81", size=(70, -1))
        self.cy_txt = wx.TextCtrl(panel, value="-1.59", size=(70, -1))
        self.cz_txt = wx.TextCtrl(panel, value="4.48", size=(70, -1))
        for label, control in (
            ("X:", self.cx_txt),
            ("Y:", self.cy_txt),
            ("Z:", self.cz_txt),
        ):
            center_box.Add(wx.StaticText(panel, label=label),
                           flag=wx.ALIGN_CENTER_VERTICAL | wx.LEFT | wx.RIGHT, border=5)
            center_box.Add(control, flag=wx.ALIGN_CENTER_VERTICAL)
        vbox.Add(center_box, flag=wx.LEFT | wx.RIGHT | wx.BOTTOM, border=10)

        # --- Box Size Input ---
        size_box = wx.BoxSizer(wx.HORIZONTAL)
        size_box.Add(wx.StaticText(panel, label="Box Size:"),
                     flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=8)
        self.sx_txt = wx.TextCtrl(panel, value="26", size=(70, -1))
        self.sy_txt = wx.TextCtrl(panel, value="26", size=(70, -1))
        self.sz_txt = wx.TextCtrl(panel, value="26", size=(70, -1))
        for label, control in (
            ("X:", self.sx_txt),
            ("Y:", self.sy_txt),
            ("Z:", self.sz_txt),
        ):
            size_box.Add(wx.StaticText(panel, label=label),
                         flag=wx.ALIGN_CENTER_VERTICAL | wx.LEFT | wx.RIGHT, border=5)
            size_box.Add(control, flag=wx.ALIGN_CENTER_VERTICAL)
        vbox.Add(size_box, flag=wx.LEFT | wx.RIGHT | wx.BOTTOM, border=10)

        # --- Convert Button ---
        self.run_btn = wx.Button(panel, label="Prepare Receptor")
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
        pdb_file = self.file_picker.GetPath().strip()

        if not pdb_file or not os.path.isfile(pdb_file):
            wx.MessageBox("Please select an existing receptor PDB file.", "Error", wx.OK | wx.ICON_ERROR)
            return

        try:
            center = [float(self.cx_txt.GetValue()),
                      float(self.cy_txt.GetValue()),
                      float(self.cz_txt.GetValue())]
            size = [float(self.sx_txt.GetValue()),
                    float(self.sy_txt.GetValue()),
                    float(self.sz_txt.GetValue())]
        except ValueError:
            wx.MessageBox("Box center and size values must be numbers.", "Error", wx.OK | wx.ICON_ERROR)
            return

        if any(v <= 0 for v in size):
            wx.MessageBox("Box size values must be greater than zero.", "Error", wx.OK | wx.ICON_ERROR)
            return

        mk_receptor = os.path.join(str(pref.AutoDock_dir), 'venv', 'bin', 'mk_prepare_receptor.py')
        if not os.path.isfile(mk_receptor):
            wx.MessageBox(
                f"mk_prepare_receptor.py not found:\n{mk_receptor}\n\nCheck AutoDock_dir in preferences.py.",
                "Error", wx.OK | wx.ICON_ERROR)
            return

        # Disable controls during execution
        self.run_btn.Disable()
        self.gauge.SetValue(0)
        self.log_txt.Clear()

        # Run conversion in a separate thread to keep the GUI active
        thread = threading.Thread(target=self.run_conversion, args=(pdb_file, mk_receptor, center, size))
        # Keep the process alive until the receptor has been fully prepared.
        thread.daemon = False
        self.worker = thread
        thread.start()

    def run_conversion(self, pdb_file, mk_receptor, center, size):
        dir_name = os.path.dirname(os.path.abspath(pdb_file))
        base_name = os.path.splitext(os.path.basename(pdb_file))[0]

        cmd = [
            mk_receptor,
            '--read_pdb', pdb_file,
            '-o', base_name,
            '-p', '-j', '-g',
            '--box_center', str(center[0]), str(center[1]), str(center[2]),
            '--box_size', str(size[0]), str(size[1]), str(size[2]),
        ]

        self.log(f"Preparing receptor: {pdb_file}")
        self.log(f"Box center: {center[0]} {center[1]} {center[2]}")
        self.log(f"Box size: {size[0]} {size[1]} {size[2]}")
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

            if process.returncode == 0:
                wx.CallAfter(self.gauge.SetValue, 100)
                self.log(f"    --> Successfully created: {os.path.join(dir_name, base_name)}.pdbqt")
                self.log("\n" + "="*50 + "\nReceptor preparation finished.")
            else:
                self.log(f"\n    --> ERROR: mk_prepare_receptor.py exited with code {process.returncode}")

        except Exception as e:
            self.log(f"    --> Failed to execute mk_prepare_receptor.py: {str(e)}")

        wx.CallAfter(self.run_btn.Enable)


if __name__ == "__main__":
    app = wx.App(False)
    frame = Pdb2PdbqtFrame()
    frame.Show()
    app.MainLoop()
