#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 17 12:10:31 2023

@author: Jan Jakubik jan.jakubik@fgu.cas.cz
"""
import os
import wx
import subprocess
import preferences as pref
install_dir = (pref.install_dir)
import pandas as pd

#error log
if os.path.isfile("temp.log"):
    os.remove("temp.log")
err = ''

#select log file
dialog = wx.FileDialog(None, "Choose log file of rescoring", os.getcwd(), "","*rescored.log", wx.FD_OPEN)
if dialog.ShowModal() == wx.ID_OK:
    log=dialog.GetPath()
    dialog.Destroy()
    #get base name
    full_path=os.path.splitext(log)[0]
    target=os.path.basename(full_path).split('/')[-1]
    target_log=(str(target)+'.log')
    target_tab=(str(target)+'.tab')
    target_hbonds=(str(target)+'_hbonds.tab')
    target_ionic=(str(target)+'_ionic.tab')
    target_cationPi=(str(target)+'_cationPi.tab')
    target_PiPi=(str(target)+'_PiPi.tab')
    target_distances=(str(target)+'_distances.tab')
    #Info
    msg = (str(target)+' poses will be sorted.')
    dialog = wx.MessageDialog(None, msg, 'Info', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()
    #backup copy of log file
    cmd0 =('cp '+str(target_log)+' '+str(target_tab))
    p0 = subprocess.Popen(cmd0, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p0.wait()
    #remove 9th line of log file
    cmd1 = ("sed -i \'9d\' "+str(target_tab))
    p1 = subprocess.Popen(cmd1, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p1.wait()
    #remove lines 1 to 7
    cmd2 = ("sed -i \'1,7d\' "+str(target_tab))
    p2 = subprocess.Popen(cmd2, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    p2.wait()
    #read in poses and energies
    poses = pd.read_csv(target_tab,delimiter='|')
    poses.drop(poses.tail(2).index,inplace=True)
    x_data = poses.loc[:,"#Lig."].astype(int)
    energy = poses.loc[:,"Bind.energy[kcal/mol]"].astype(float)
    table = pd.concat([x_data,energy], ignore_index=False, axis=1)
    #H-bonds
    if os.path.isfile(target_hbonds):
        hbonds = pd.read_csv(target_hbonds,delimiter='|')
        y_data = hbonds.loc[:," #____Pose"].astype(int)
        z_data = hbonds.loc[:," HBond Residues"]
        yz_data = pd.concat([y_data,z_data], ignore_index=False, axis=1)
        table = pd.merge(table, yz_data, left_on='#Lig.',right_on=' #____Pose')
        table = table.drop([' #____Pose'], axis=1)
    #inonic bonds
    if os.path.isfile(target_ionic):
        ionic = pd.read_csv(target_ionic,delimiter='|')
        y_data = ionic.loc[:," #____Pose"].astype(int)
        z_data = ionic.loc[:," IBond Residues"]
        yz_data = pd.concat([y_data,z_data], ignore_index=False, axis=1)
        table = pd.merge(table, yz_data, left_on='#Lig.',right_on=' #____Pose')
        table = table.drop([' #____Pose'], axis=1)
    #cation-Pi
    if os.path.isfile(target_cationPi):
        cationPi = pd.read_csv(target_cationPi,delimiter='|')
        y_data = cationPi.loc[:," #____Pose"].astype(int)
        z_data = cationPi.loc[:," CatPi Residues"]
        yz_data = pd.concat([y_data,z_data], ignore_index=False, axis=1)
        table = pd.merge(table, yz_data, left_on='#Lig.',right_on=' #____Pose')
        table = table.drop([' #____Pose'], axis=1)
    #Pi-Pi
    if os.path.isfile(target_PiPi):
        PiPi = pd.read_csv(target_PiPi,delimiter='|')
        y_data = PiPi.loc[:," #____Pose"].astype(int)
        z_data = PiPi.loc[:," PiPi Residues"]
        yz_data = pd.concat([y_data,z_data], ignore_index=False, axis=1)
        table = pd.merge(table, yz_data, left_on='#Lig.',right_on=' #____Pose')
        table = table.drop([' #____Pose'], axis=1)
    #distances
    if os.path.isfile(target_distances):
        distances = pd.read_csv(target_distances,delimiter='|')
        y_data = distances.loc[:," #____Pose"].astype(int)
        z_data = distances.iloc[:,1]
        yz_data = pd.concat([y_data,z_data], ignore_index=False, axis=1)
        table = pd.merge(table, yz_data, left_on='#Lig.',right_on=' #____Pose')
        table = table.drop([' #____Pose'], axis=1)
    #write output
    output=(str(target)+'_sorted.csv')
    table.to_csv(output)
    #End
    dialog = wx.MessageDialog(None, 'Sorting poses finished.', 'Done', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()
    err = "Done"
else:
    msg = ('No log file selected!\nExiting!')
    err = "Canceled"
    dialog = wx.MessageDialog(None, msg, 'Error', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()

log = open("temp.log","w")
log.write(err)
log.close()
