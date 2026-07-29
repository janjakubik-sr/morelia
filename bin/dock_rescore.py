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

#select sce file
dialog = wx.FileDialog(None, "Choose input sce file", os.getcwd(), "","*.sce", wx.FD_OPEN)
if dialog.ShowModal() == wx.ID_OK:
    sce=dialog.GetPath()
    #get base name
    full_path=os.path.splitext(sce)[0]
    target=os.path.basename(full_path).split('/')[-1]
    msg = ('Rescoring will be performed on\n\n'+str(target))
    dialog = wx.MessageDialog(None, msg, 'Info', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()
    #commands
    cmd0 = ('cp '+str(install_dir)+'/YASARA/dock_rescore.mcr ./')
    p0 = subprocess.Popen(cmd0, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p0.wait()
    cmd1 = ("sed -i \'s/macro_target/"+str(target)+"/g\' dock_rescore.mcr")
    p1 = subprocess.Popen(cmd1, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p1.wait()
    #Run now?
    msg = ('Rescoring macro written to dock_rescore.mcr\nRun it now?')
    dialog = wx.MessageDialog(None, msg, 'Info', wx.YES_NO | wx.CANCEL | wx.NO_DEFAULT)
    if dialog.ShowModal() == wx.ID_YES:
        cmd6 = ('yasara -txt dock_rescore.mcr > rescore.log &')
        p6 = subprocess.Popen(cmd6, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p6.wait()
        msg = ('Rescoring of binding poses was started in the background.')
        message = wx.MessageDialog(None, msg, 'Done', wx.OK)
        message.ShowModal()
        message.Destroy()
        err = msg
    else:
        msg = ('Rescoring can be run from comman line.')
        message = wx.MessageDialog(None, msg, 'Done', wx.OK)
        message.ShowModal()
        message.Destroy()
        err = msg
    dialog.Destroy()
else:
    msg = ('No sce input file selected!\nExiting!')
    err = "Canceled"
    dialog = wx.MessageDialog(None, msg, 'Error', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()

log = open("temp.log","w")
log.write(err)
log.close()    
