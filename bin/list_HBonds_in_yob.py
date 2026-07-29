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
dialog = wx.FileDialog(None, "Choose input yob file", os.getcwd(), "","*.yob", wx.FD_OPEN)
if dialog.ShowModal() == wx.ID_OK:
    yob=dialog.GetPath()
    #get base name
    full_path=os.path.splitext(yob)[0]
    target=os.path.basename(full_path).split('/')[-1]
    msg = ('HBonds in\n'+str(target)+'\nwill be listed')
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
    #commands
    cmd0 = ('cp '+str(install_dir)+'/YASARA/list_HBonds_in_yob.mcr ./')
    p0 = subprocess.Popen(cmd0, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p0.wait()
    cmd1 = ("sed -i \'s/yob_file/"+str(target)+"/g\' list_HBonds_in_yob.mcr")
    p1 = subprocess.Popen(cmd1, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p1.wait()
    cmd3 = ("sed -i \'s/mol_receptor/"+str(rec)+"/g\' list_HBonds_in_yob.mcr")
    p3 = subprocess.Popen(cmd3, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p3.wait()
    cmd4 = ("sed -i \'s/mol_ligand/"+str(lig)+"/g\' list_HBonds_in_yob.mcr")
    p4 = subprocess.Popen(cmd4, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p4.wait()
    cmd5 = ('yasara -txt list_HBonds_in_yob.mcr')
    p5 = subprocess.Popen(cmd5, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p5.wait()
    msg = ('Ligand-receptor HBonds in\n'+str(target)+'\nare listed in console')
    dialog = wx.MessageDialog(None, msg, 'Info', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()
else:
    msg = ('No yob input file selected!\nExiting!')
    err = "Canceled"
    dialog = wx.MessageDialog(None, msg, 'Error', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()

log = open("temp.log","w")
log.write(err)
log.close()
