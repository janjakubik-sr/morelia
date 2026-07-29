#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spyder Editor

Archive current folder.
"""

#No edits below

import glob
import os
import wx
import subprocess
import preferences as pref
install_dir = (pref.install_dir)

#error log
if os.path.isfile("temp.log"):
    os.remove("temp.log")
err = ''

#select sce file
dialog = wx.FileDialog(None, "Choose base name file", os.getcwd(), "","*.sce", wx.FD_OPEN)
if dialog.ShowModal() == wx.ID_OK:
    sce=dialog.GetPath()
    #get base name
    full_path=os.path.splitext(sce)[0]
    target=os.path.basename(full_path).split('/')[-1]
    msg = (str(target)+'* files will compressed to archive')
    dialog = wx.MessageDialog(None, msg, 'Info', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()
    #commands
    cmd0 = ('tar -cJvf '+str(target)+'.tar.xz '+str(target)+'*')
    p0 = subprocess.Popen(cmd0, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    dialog = wx.MessageDialog(None, 'Compressing files\nPlease wait', 'Please wait', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()
    p0.wait()
    cmd1 = ('rm '+str(target)+'*.yob')
    p1 = subprocess.Popen(cmd1, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p1.wait()
    cmd2 = ('rm '+str(target)+'*.adr')
    p2 = subprocess.Popen(cmd2, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p2.wait()
    cmd3 = ('rm '+str(target)+'*.map')
    p3 = subprocess.Popen(cmd3, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p3.wait()
    cmd4 = ('rm '+str(target)+'*.fld')
    p4 = subprocess.Popen(cmd4, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p4.wait()
    cmd5 = ('rm rescored.log')
    p5 = subprocess.Popen(cmd5, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p5.wait()
    cmd6 = ('rm *check*')
    p6 = subprocess.Popen(cmd6, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p6.wait()
    #Finished
    dialog = wx.MessageDialog(None, 'Archiving finnished.', 'Done', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()
    err = "Done."
else:
    msg = ('Archiving canceled!\nExiting!')
    err = "Canceled"
    dialog = wx.MessageDialog(None, msg, 'Error', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()

log = open("temp.log","w")
log.write(err)
log.close()
