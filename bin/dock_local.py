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
    file_name=os.path.basename(full_path).split('/')[-1]
    target = file_name.replace("_complex","")
    msg = ('Local search will be performed on\n\n'+str(target))
    dialog = wx.MessageDialog(None, msg, 'Info', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()
    #number of runs
    ask = ('Give number of runs')
    dlg = wx.TextEntryDialog(None, ask,"Number of runs","888")
    if dlg.ShowModal() == wx.ID_OK:
        nruns = dlg.GetValue()
    dlg.Destroy()
    #RMSD min
    ask = ('Give minimum RMSD')
    dlg = wx.TextEntryDialog(None, ask,"RMSD min","1.0")
    if dlg.ShowModal() == wx.ID_OK:
        rmsd = dlg.GetValue()
    dlg.Destroy()
    #commands
    cmd0 = ('cp '+str(install_dir)+'/YASARA/dock_runlocal.mcr ./')
    p0 = subprocess.Popen(cmd0, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p0.wait()
    cmd1 = ("sed -i \'s/macro_target/"+str(target)+"/g\' dock_runlocal.mcr")
    p1 = subprocess.Popen(cmd1, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p1.wait()
    cmd2 = ("sed -i \'s/number_of_runs/"+str(nruns)+"/g\' dock_runlocal.mcr")
    p2 = subprocess.Popen(cmd2, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p2.wait()
    cmd3 = ("sed -i \'s/rmsd_min/"+str(rmsd)+"/g\' dock_runlocal.mcr")
    p3 = subprocess.Popen(cmd3, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p3.wait()
    cmd6 = ('yasara -txt dock_runlocal.mcr')
    p6 = subprocess.Popen(cmd6, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    #Wait
    dialog = wx.MessageDialog(None, 'Local search docking\nPlease wait', 'Please wait', wx.OK)
    dialog.ShowModal()
    p6.wait()
    dialog.Destroy()
    #Docking finished
    dialog = wx.MessageDialog(None, 'Local search docking finnished.', 'Done', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()
    err ="Done"
else:
    msg = ('No sce input file selected!\nExiting!')
    err = "Canceled"
    dialog = wx.MessageDialog(None, msg, 'Error', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()

log = open("temp.log","w")
log.write(err)
log.close()