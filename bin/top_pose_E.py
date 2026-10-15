#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 09:45:10 2024

@author: roshi
"""
import glob, os
import pandas as pd
import wx

#error log
if os.path.isfile("temp.log"):
    os.remove("temp.log")
err = ''

table = open("top_poses_E.csv","w")

for file in sorted(glob.glob('*/*_sorted.csv')):
    basename=os.path.basename(file)
    base=os.path.splitext(basename)[0]
    poses = pd.read_csv(file,delimiter=',')
    top_pose_E = poses.loc[0,"Bind.energy[kcal/mol]"].astype(float)
    table.write(str(base)+','+str(top_pose_E)+'\n')

table.close()
dialog = wx.MessageDialog(None, 'Extracting energies finished.', 'Done', wx.OK)
dialog.ShowModal()
dialog.Destroy()
err = "Done"

log = open("temp.log","w")
log.write(err)
log.close()
