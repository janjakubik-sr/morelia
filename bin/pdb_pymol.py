#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 14:15:28 2023

@author: Jan Jakubik jan.jakubik@fgu.cas.cz
"""
import os
import wx
import subprocess

#error log
if os.path.isfile("temp.log"):
    os.remove("temp.log")
err = ''

#select pdf file
dialog = wx.FileDialog(None, "Choose yob file to convert", os.getcwd(), "","*.pdb", wx.FD_OPEN)
if dialog.ShowModal() == wx.ID_OK:
    pdb=dialog.GetPath()
    dialog.Destroy()
    cmd0 = ("pymol -x "+str(pdb))
    p0 = subprocess.Popen(cmd0, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
else:
    msg = ('No file selected!\nExiting!')
    err = msg
    dialog = wx.MessageDialog(None, msg, 'Error', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()

log = open("temp.log","w")
log.write(err)
log.close()