#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 11:58:47 2026

@author: roshi
"""
import os
import glob
import csv
import subprocess
import threading
import wx

class ExtractPosesFrame(wx.Frame):
    def __init__(self, parent, title="Top Poses Extractor"):
        super().__init__(parent, title=title, size=(600, 480))
        self.InitUI()
        self.Centre()

    def InitUI(self):
        panel = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        # Directory Glob Input
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        lbl_glob = wx.StaticText(panel, label="Directory Glob:")
        self.txt_glob = wx.TextCtrl(panel, value="M5_*")
        hbox1.Add(lbl_glob, flag=wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, border=8)
        hbox1.Add(self.txt_glob, proportion=1)
        vbox.Add(hbox1, flag=wx.EXPAND | wx.ALL, border=10)

        # Action Button
        self.btn_run = wx.Button(panel, label="Run Extraction")
        self.btn_run.Bind(wx.EVT_BUTTON, self.OnRun)
        vbox.Add(self.btn_run, flag=wx.ALIGN_RIGHT | wx.RIGHT | wx.BOTTOM, border=10)

        # Output Log Text Area
        lbl_log = wx.StaticText(panel, label="Log Output:")
        vbox.Add(lbl_log, flag=wx.LEFT | wx.BOTTOM, border=5)
        self.txt_log = wx.TextCtrl(panel, style=wx.TE_MULTILINE | wx.TE_READONLY)
        vbox.Add(self.txt_log, proportion=1, flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, border=10)

        panel.SetSizer(vbox)

    def Log(self, msg):
        wx.CallAfter(self.txt_log.AppendText, msg + "\n")

    def OnRun(self, event):
        glob_pattern = self.txt_glob.GetValue().strip()
        if not glob_pattern:
            wx.MessageBox("Please enter a valid directory pattern.", "Error", wx.OK | wx.ICON_ERROR)
            return

        self.btn_run.Disable()
        self.txt_log.Clear()
        
        # Run processing in a background thread to keep GUI responsive
        threading.Thread(target=self.ProcessDirectories, args=(glob_pattern,), daemon=True).start()

    def ProcessDirectories(self, pattern):
        matching_dirs = [d for d in glob.glob(pattern) if os.path.isdir(d)]

        if not matching_dirs:
            self.Log(f"No directories found matching pattern: '{pattern}'")
            wx.CallAfter(self.btn_run.Enable)
            return

        self.Log(f"Found {len(matching_dirs)} directory/directories to process...\n")

        for d in matching_dirs:
            folder_name = os.path.basename(os.path.normpath(d))
            
            # Find any *_rescored_sorted.csv file in the folder
            csv_matches = glob.glob(os.path.join(d, "*rescored_sorted.csv"))
            if not csv_matches:
                self.Log(f"[SKIP] No '*_rescored_sorted.csv' found in: {d}")
                continue
            csv_path = csv_matches[0]

            # Find any *.tar.xz archive in the folder
            tar_matches = glob.glob(os.path.join(d, "*.tar.xz"))
            if not tar_matches:
                self.Log(f"[SKIP] No '*.tar.xz' archive found in: {d}")
                continue
            tar_path = tar_matches[0]
            tar_filename = os.path.basename(tar_path)

            # Extract rank value from 2nd line, 2nd column
            try:
                with open(csv_path, 'r', newline='') as f:
                    reader = csv.reader(f)
                    _ = next(reader, None)  # Skip header line
                    second_row = next(reader, None)

                if not second_row or len(second_row) < 2:
                    self.Log(f"[ERROR] {csv_path} does not have a 2nd row with at least 2 columns.")
                    continue

                raw_val = second_row[1].strip()
                val = int(raw_val)
                pose_id = f"{val:03d}"
            except Exception as e:
                self.Log(f"[ERROR] Failed to read pose ID from {csv_path}: {e}")
                continue

            # Target pattern for tar extraction (matches pose ID suffix)
            target_pattern = f"*rescored*{pose_id}.yob"

            # Execute extraction
            cmd = ["tar", "-xvf", tar_filename, "--wildcards", target_pattern]
            self.Log(f"[PROCESSING] {folder_name}: Extracting target '{target_pattern}' from '{tar_filename}'...")

            try:
                res = subprocess.run(cmd, cwd=d, capture_output=True, text=True)
                if res.returncode == 0:
                    out = res.stdout.strip()
                    self.Log(out if out else "[SUCCESS] Extracted successfully.")
                else:
                    self.Log(f"[ERROR] tar command failed in {d}:\n{res.stderr.strip()}")
            except Exception as e:
                self.Log(f"[ERROR] Execution failed for {folder_name}: {e}")

        self.Log("\n--- Processing Complete ---")
        wx.CallAfter(self.btn_run.Enable)

if __name__ == '__main__':
    app = wx.App(False)
    frame = ExtractPosesFrame(None)
    frame.Show()
    app.MainLoop()