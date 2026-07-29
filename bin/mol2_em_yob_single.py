#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 10:13:16 2023

@author: Jan Jakubik jan.jakubik@fgu.cas.cz
"""
import os
import wx
import subprocess
import preferences as pref
install_dir = (pref.install_dir)


#select YOb file
dialog = wx.FileDialog(None, "Choose input mol2 file", os.getcwd(), "","*.mol2", wx.FD_OPEN)
if dialog.ShowModal() == wx.ID_OK:
    path=dialog.GetPath()
dialog.Destroy()

cmd0 = ('cp '+str(install_dir)+'/YASARA/mol2_em_yob.mcr temp.mcr')
p0 = subprocess.Popen(cmd0, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
p0.wait()
file=os.path.basename(path)
base=os.path.splitext(file)[0]
cmd1 = ("sed -i \'s/target/"+str(file)+"/g\' temp.mcr")
p1 = subprocess.Popen(cmd1, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
p1.wait()
cmd2 = ("sed -i \'s/MyName/L/g\' temp.mcr")
p2 = subprocess.Popen(cmd2, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
p2.wait()
output=(str(base)+'.yob')
cmd3 = ("sed -i \'s/output/"+str(output)+"/g\' temp.mcr")
p3 = subprocess.Popen(cmd3, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
p3.wait()
cmd4 = ("yasara -txt temp.mcr")
p4 = subprocess.Popen(cmd4, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
p4.wait()