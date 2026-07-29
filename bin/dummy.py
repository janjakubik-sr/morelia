#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 12:35:07 2016

@author: Jan Jakubik jan.jakubik@fgu.cas.cz
"""
import glob,os
import subprocess
import wx
import preferences as pref
install_dir = (pref.install_dir)

#error log
if os.path.isfile("temp.log"):
    os.remove("temp.log")
msg = ''

#test wheter we have files
count = len(glob.glob('desmond_md_job_*/desmond_md_job_*.sh'))
if count == 0:
    msg = ('Aborted!\nNo MD job files found!')
    dialog = wx.MessageDialog(None, msg, 'Error', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()
else:
    for file in sorted(glob.glob('*.dat')):
        # run dummy
        cmd = (str(install_dir)+'dummy.sh')
        p1 = subprocess.Popen(cmd, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p1.wait()
        msg = ('MD job files edited')
        dialog = wx.MessageDialog(None, msg, 'Info', wx.OK)
        dialog.ShowModal()
        dialog.Destroy()

log = open("temp.log","w")
log.write(msg)
log.close()