#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 14:38:15 2023

@author: Jan Jakubik jan.jakubik@fgu.cas.cz
"""

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
    msg = ('Ligand-receptor intaractions in\n\n'+str(target)+'\n\nwill be analysed.')
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

#Select interactions to analyze
lst = [ 'Hydrogen bonds', 'Hydrophobic interactions', 'Ionic bonds', 'Pi-cation', 'PiPi']
dlg = wx.MultiChoiceDialog(None, "Select interaction types to analyze.", "Select", lst)
if (dlg.ShowModal() == wx.ID_OK):
    selections = dlg.GetSelections()
    dlg.Destroy()
    if 0 in selections:
        #Hbonds
        cmd0 = ('cp '+str(install_dir)+'/YASARA/dock_hbonds.mcr ./')
        p0 = subprocess.Popen(cmd0, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p0.wait()
        cmd1 = ("sed -i \'s/macro_target/"+str(target)+"/g\' dock_hbonds.mcr")
        p1 = subprocess.Popen(cmd1, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p1.wait()
        cmd2 = ("sed -i \'s/mol_receptor/"+str(rec)+"/g\' dock_hbonds.mcr")
        p2 = subprocess.Popen(cmd2, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p2.wait()
        cmd3 = ("sed -i \'s/mol_ligand/"+str(lig)+"/g\' dock_hbonds.mcr")
        p3 = subprocess.Popen(cmd3, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p3.wait()
        cmd4 = ("sed -i \'s/first_pose/"+str(fp)+"/g\' dock_hbonds.mcr")
        p4 = subprocess.Popen(cmd4, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p4.wait()
        cmd5 = ("sed -i \'s/last_pose/"+str(lp)+"/g\' dock_hbonds.mcr")
        p5 = subprocess.Popen(cmd5, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p5.wait()
        cmd6 = ('yasara -txt dock_hbonds.mcr')
        p6 = subprocess.Popen(cmd6, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        #Wait
        dialog = wx.MessageDialog(None, 'Analysing hydrogen bonds\nPlease wait', 'Please wait', wx.OK)
        dialog.ShowModal()
        p6.wait()
        dialog.Destroy()
    if 1 in selections:
        #hydrophobic
        cmd0 = ('cp '+str(install_dir)+'/YASARA/dock_hydrophobic.mcr ./')
        p0 = subprocess.Popen(cmd0, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p0.wait()
        cmd1 = ("sed -i \'s/macro_target/"+str(target)+"/g\' dock_hydrophobic.mcr")
        p1 = subprocess.Popen(cmd1, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p1.wait()
        cmd2 = ("sed -i \'s/mol_receptor/"+str(rec)+"/g\' dock_hydrophobic.mcr")
        p2 = subprocess.Popen(cmd2, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p2.wait()
        cmd3 = ("sed -i \'s/mol_ligand/"+str(lig)+"/g\' dock_hydrophobic.mcr")
        p3 = subprocess.Popen(cmd3, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p3.wait()
        cmd4 = ("sed -i \'s/first_pose/"+str(fp)+"/g\' dock_hydrophobic.mcr")
        p4 = subprocess.Popen(cmd4, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p4.wait()
        cmd5 = ("sed -i \'s/last_pose/"+str(lp)+"/g\' dock_hydrophobic.mcr")
        p5 = subprocess.Popen(cmd5, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p5.wait()
        cmd6 = ('yasara -txt dock_hydrophobic.mcr')
        p6 = subprocess.Popen(cmd6, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        #Wait
        dialog = wx.MessageDialog(None, 'Analysing hyrophobic interactions\nPlease wait', 'Please wait', wx.OK)
        dialog.ShowModal()
        p6.wait()
    if 2 in selections:
        #ionic
        cmd0 = ('cp '+str(install_dir)+'/YASARA/dock_ionic.mcr ./')
        p0 = subprocess.Popen(cmd0, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p0.wait()
        cmd1 = ("sed -i \'s/macro_target/"+str(target)+"/g\' dock_ionic.mcr")
        p1 = subprocess.Popen(cmd1, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p1.wait()
        cmd2 = ("sed -i \'s/mol_receptor/"+str(rec)+"/g\' dock_ionic.mcr")
        p2 = subprocess.Popen(cmd2, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p2.wait()
        cmd3 = ("sed -i \'s/mol_ligand/"+str(lig)+"/g\' dock_ionic.mcr")
        p3 = subprocess.Popen(cmd3, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p3.wait()
        cmd4 = ("sed -i \'s/first_pose/"+str(fp)+"/g\' dock_ionic.mcr")
        p4 = subprocess.Popen(cmd4, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p4.wait()
        cmd5 = ("sed -i \'s/last_pose/"+str(lp)+"/g\' dock_ionic.mcr")
        p5 = subprocess.Popen(cmd5, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p5.wait()
        cmd6 = ('yasara -txt dock_ionic.mcr')
        p6 = subprocess.Popen(cmd6, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        #Wait
        dialog = wx.MessageDialog(None, 'Analysing ionic interactions\nPlease wait', 'Please wait', wx.OK)
        dialog.ShowModal()
        p6.wait()
        dialog.Destroy()
    if 3 in selections:
        #cation-Pi
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
    if 4 in selections:
        #Pi-Pi
        cmd0 = ('cp '+str(install_dir)+'/YASARA/dock_PiPi.mcr ./')
        p0 = subprocess.Popen(cmd0, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p0.wait()
        cmd1 = ("sed -i \'s/macro_target/"+str(target)+"/g\' dock_PiPi.mcr")
        p1 = subprocess.Popen(cmd1, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p1.wait()
        cmd2 = ("sed -i \'s/mol_receptor/"+str(rec)+"/g\' dock_PiPi.mcr")
        p2 = subprocess.Popen(cmd2, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p2.wait()
        cmd3 = ("sed -i \'s/mol_ligand/"+str(lig)+"/g\' dock_PiPi.mcr")
        p3 = subprocess.Popen(cmd3, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p3.wait()
        cmd4 = ("sed -i \'s/first_pose/"+str(fp)+"/g\' dock_PiPi.mcr")
        p4 = subprocess.Popen(cmd4, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p4.wait()
        cmd5 = ("sed -i \'s/last_pose/"+str(lp)+"/g\' dock_PiPi.mcr")
        p5 = subprocess.Popen(cmd5, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        p5.wait()
        cmd6 = ('yasara -txt dock_PiPi.mcr')
        p6 = subprocess.Popen(cmd6, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        #Wait
        dialog = wx.MessageDialog(None, 'Analysing PiPi interactions\nPlease wait', 'Please wait', wx.OK)
        dialog.ShowModal()
        p6.wait()
        dialog.Destroy()
    #Analysis finished
    dialog = wx.MessageDialog(None, 'All analyses finnished.', 'Done', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()
    strings = [lst[x] for x in selections]
    log = open("temp.log","w")
    log.write("Intaractions analyzed:\n")
    for j in strings:
        log.write(str(j)+'\n')
    log.close()
else:
    msg = ('Analysis canceled!\nExiting!')
    err = "Canceled"
    dialog = wx.MessageDialog(None, msg, 'Error', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()
    log = open("temp.log","w")
    log.write(err)
    log.close()
