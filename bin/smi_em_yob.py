#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 10:13:16 2023

@author: Jan Jakubik jan.jakubik@fgu.cas.cz
"""
import glob, os
import wx
import subprocess
import preferences as pref
install_dir = (pref.install_dir)
i = 0

count = len(glob.glob('*.smi'))
if count == 0:
    msg = ('No .smi files found!\n')
    err = msg
    dialog = wx.MessageDialog(None, msg, 'Error', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()
else:
    for file in sorted(glob.glob('*.smi')):
        i = i + 1
        cmd0 = ('cp '+str(install_dir)+'/YASARA/smi_em_yob.mcr temp.mcr')
        p0 = subprocess.Popen(cmd0, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p0.wait()
        cmd1 = ("sed -i \'s/target/"+str(file)+"/g\' temp.mcr")
        p1 = subprocess.Popen(cmd1, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p1.wait()
        cmd2 = ("sed -i \'s/MyName/L"+str(i)+"/g\' temp.mcr")
        p2 = subprocess.Popen(cmd2, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p2.wait()
        base=os.path.splitext(file)[0]
        output=(str(base)+'.yob')
        cmd3 = ("sed -i \'s/output/"+str(output)+"/g\' temp.mcr")
        p3 = subprocess.Popen(cmd3, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p3.wait()
        cmd4 = ("yasara -txt temp.mcr")
        p4 = subprocess.Popen(cmd4, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p4.wait()
