# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
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
dialog = wx.FileDialog(None, "Choose input sce file", os.getcwd(), "","*.sce", wx.FD_OPEN)
if dialog.ShowModal() == wx.ID_OK:
    sce=dialog.GetPath()
    #get base name
    full_path=os.path.splitext(sce)[0]
    target=os.path.basename(full_path).split('/')[-1]
    msg = ('Cation-Pi interaction analysis will be performed on\n\n'+str(target))
    dialog = wx.MessageDialog(None, msg, 'Info', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()
    #receptor
    ask = ('Give ASL for receptor e.g. Mol A')
    dlg = wx.TextEntryDialog(None, ask,"Receptor","Mol A")
    if dlg.ShowModal() == wx.ID_OK:
        rec = dlg.GetValue()
    dlg.Destroy()
    #ligand
    ask = ('Give ASL for ligand e.g. Mol L')
    dlg = wx.TextEntryDialog(None, ask,"Ligand","Mol L")
    if dlg.ShowModal() == wx.ID_OK:
        lig = dlg.GetValue()
    dlg.Destroy()
    #poses
    count = len(glob.glob(str(target)+'_???.yob'))
    fp = 1
    lp = count
    #commands
    cmd0 = ('cp '+str(install_dir)+'/YASARA/dock_cationPi.mcr ./')
    p0 = subprocess.Popen(cmd0, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p0.wait()
    cmd1 = ("sed -i \'s/macro_target/"+str(target)+"/g\' dock_cationPi.mcr")
    p1 = subprocess.Popen(cmd1, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p1.wait()
    cmd2 = ("sed -i \'s/mol_receptor/"+str(rec)+"/g\' dock_cationPi.mcr")
    p2 = subprocess.Popen(cmd2, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p2.wait()
    cmd3 = ("sed -i \'s/mol_ligand/"+str(lig)+"/g\' dock_cationPi.mcr")
    p3 = subprocess.Popen(cmd3, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p3.wait()
    cmd4 = ("sed -i \'s/first_pose/"+str(fp)+"/g\' dock_cationPi.mcr")
    p4 = subprocess.Popen(cmd4, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p4.wait()
    cmd5 = ("sed -i \'s/last_pose/"+str(lp)+"/g\' dock_cationPi.mcr")
    p5 = subprocess.Popen(cmd5, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p5.wait()
    cmd6 = ('yasara -txt dock_cationPi.mcr')
    p6 = subprocess.Popen(cmd6, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    #Wait
    dialog = wx.MessageDialog(None, 'Analysing caton-Pi interactions\nPlease wait', 'Please wait', wx.OK)
    dialog.ShowModal()
    p6.wait()
    dialog.Destroy()
    #Analysis finished
    dialog = wx.MessageDialog(None, 'Analysis of cation-Pi interactions finnished.', 'Done', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()
    err = "Done"
else:
    msg = ('No sce input file selected!\nExiting!')
    err = "Canceled"
    dialog = wx.MessageDialog(None, msg, 'Error', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()

log = open("temp.log","w")
log.write(err)
log.close()