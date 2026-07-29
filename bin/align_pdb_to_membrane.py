#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#No edits below

import os
import wx
import subprocess
import preferences as pref
install_dir = (pref.install_dir)

#error log
if os.path.isfile("temp.log"):
    os.remove("temp.log")
err = ''

#select yob file
dialog = wx.FileDialog(None, "Choose pdb file to align", os.getcwd(), "","*.pdb", wx.FD_OPEN)
if dialog.ShowModal() == wx.ID_OK:
    yob=dialog.GetPath()
else:
    msg = ('No file selected!\nExiting!')
    err = msg
    dialog = wx.MessageDialog(None, msg, 'Error', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()

#get base name
full_path=os.path.splitext(yob)[0]
target=os.path.basename(full_path).split('/')[-1]

msg = (str(target)+'.pdb\nwill be aligned to membrane')
dialog = wx.MessageDialog(None, msg, 'Info', wx.OK)
dialog.ShowModal()
dialog.Destroy()

cmd0 = ('cp '+str(install_dir)+'/YASARA/align_pdb_to_membrane.mcr ./')
p0 = subprocess.Popen(cmd0, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
p0.wait()

cmd1 = ('cp '+str(install_dir)+'/PDB/M2_3UON_OPM.pdb ./')
p1 = subprocess.Popen(cmd1, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
p1.wait()

cmd2 = ("sed -i \'s/macro_target/"+str(target)+"/g\' align_pdb_to_membrane.mcr")
p2 = subprocess.Popen(cmd2, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
p2.wait()

cmd6 = ('yasara -txt align_pdb_to_membrane.mcr')
p6 = subprocess.Popen(cmd6, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
p6.wait()

cmd7 = ('rm ./M2_3UON_OPM.pdb')
p7 = subprocess.Popen(cmd7, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
p7.wait()

#End
dialog = wx.MessageDialog(None, 'PDB aligned to membrane.', 'Done', wx.OK)
dialog.ShowModal()
dialog.Destroy()

log = open("temp.log","w")
log.write(err)
log.close()
