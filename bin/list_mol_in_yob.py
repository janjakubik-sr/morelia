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
dialog = wx.FileDialog(None, "Choose input yob file", os.getcwd(), "","*.yob", wx.FD_OPEN)
if dialog.ShowModal() == wx.ID_OK:
    yob=dialog.GetPath()
    #get base name
    full_path=os.path.splitext(yob)[0]
    target=os.path.basename(full_path).split('/')[-1]
    msg = ('Molecules in\n'+str(target)+'\nwill be listed')
    dialog = wx.MessageDialog(None, msg, 'Info', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()
    #commands
    cmd0 = ('cp '+str(install_dir)+'/YASARA/list_mol_in_yob.mcr ./')
    p0 = subprocess.Popen(cmd0, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p0.wait()
    cmd1 = ("sed -i \'s/yob_file/"+str(target)+"/g\' list_mol_in_yob.mcr")
    p1 = subprocess.Popen(cmd1, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p1.wait()
    cmd2 = ('yasara -txt list_mol_in_yob.mcr > mol_list.res')
    p2 = subprocess.Popen(cmd2, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p2.wait()
    msg = ('Molecules in\n'+str(target)+'\nare listed in console')
    dialog = wx.MessageDialog(None, msg, 'Info', wx.OK)
    dialog.ShowModal()
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