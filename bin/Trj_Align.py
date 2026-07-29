# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#No edits below

import os
import wx
import subprocess

#select cms file
dialog = wx.FileDialog(None, "Choose input cms file", os.getcwd(), "","*.cms", wx.FD_OPEN)
if dialog.ShowModal() == wx.ID_OK:
    cms=dialog.GetPath()
dialog.Destroy()

#select trajectory
dlg = wx.DirDialog(None, "Choose trajectory directory:",
                   style=wx.DD_DEFAULT_STYLE
                   | wx.DD_DIR_MUST_EXIST
                   )
if dlg.ShowModal() == wx.ID_OK:
    trj = dlg.GetPath()
dlg.Destroy()
    
#Enter output name
ask = ('Give the output base name e.g. aligned')
dlg = wx.TextEntryDialog(None, ask,"Output","aligned")
if dlg.ShowModal() == wx.ID_OK:
    out = dlg.GetValue()
dlg.Destroy()

#Enter ASL
asl = ('Give ASL for alignment e.g. chain P')
dlg = wx.TextEntryDialog(None, ask,"Output","chain P")
if dlg.ShowModal() == wx.ID_OK:
    asl = dlg.GetValue()
dlg.Destroy()

cmd = ('$SCHRODINGER/run trj_align.py -asl \'('+str(asl)+')\' '+str(cms)+' '+str(trj)+' '+str(out))
p1 = subprocess.Popen(cmd, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)

#Wait
dialog = wx.MessageDialog(None, 'Aligning trajectory\nPlease wait', 'Please wait', wx.OK)
dialog.ShowModal()
p1.wait()
dialog.Destroy()

#Align finished
dialog = wx.MessageDialog(None, 'Trajectory alignment finnished.', 'Done', wx.OK)
dialog.ShowModal()
dialog.Destroy()