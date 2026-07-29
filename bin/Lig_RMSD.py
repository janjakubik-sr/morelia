#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spyder Editor

"""

#No edits below

import os
import wx
import numpy as np

#error log
if os.path.isfile('temp.log'):
    os.remove('temp.log')
err = ''

#Extract data
frame = []
rmsd = []

if os.path.isfile('PL_RMSD.dat'):
    with open('PL_RMSD.dat') as LIG:
        for line in LIG:
            fields = line.strip().split()
            frame.append(fields[0])
            rmsd.append(fields[5])
else:
    msg = ('PL_RMSD.dat not found!\nExiting!')
    err = msg
    dialog = wx.MessageDialog(None, msg, 'Error', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()
    

lig_rmsd = np.array([0,0])
for i in range(1, 1000):
    n_frame = frame[i]
    n_rmsd = rmsd[i]
    exec('lig_rmsd_'+str(i)+'= np.stack((n_frame, n_rmsd), axis=-1)')
    exec('lig_rmsd = np.vstack((lig_rmsd,lig_rmsd_'+str(i)+'))')

#Save data file
lig_rmsd = np.delete(lig_rmsd,0,0)

#Convert to numbers
ask = ('Enter timestep in ns\n')
dlg = wx.TextEntryDialog(None, ask,"Input","0.12")
if dlg.ShowModal() == wx.ID_OK:
    ts = dlg.GetValue()
dlg.Destroy()

ts = float(ts)
n_rows = len(lig_rmsd[:,0])

for j in range(0,n_rows):
    lig_rmsd[j,0] = float(lig_rmsd[j,0])*ts

with open('Lig_RMSD.dat','w') as res:
    np.savetxt(res,lig_rmsd,fmt='%s')

data = np.loadtxt('Lig_RMSD.dat')

#Plot with matplotlib
import matplotlib.pyplot as plt
x_max = np.amax(data[:,0]) * 1.1
y_max = np.amax(data[:,1]) * 1.1
plt.title("Ligand RMSD")
plt.axis([0,x_max,0,y_max])
plt.plot(data[:,0],data[:,1])
plt.xlabel('Time [ns]')
plt.ylabel('RMSD [A]')
plt.show()

#Plot with Grace
#prepare file for Grace plot
if os.path.isfile("Lig_RMSD.draw"):
    os.remove("Lig_RMSD.draw")
draw = open("Lig_RMSD.draw","w")

ax_min = 0
ax_max = x_max
ay_min = 0
ay_max = y_max

draw.write("read xy \"Lig_RMSD.dat\" \n")
draw.write("s0 symbol 0 \n")
draw.write("s0 line type 1 \n")
draw.write("s0 line color 1 \n")

draw.write("page layout portrait \n")
draw.write("view xmin 0.2 \n")
draw.write("view xmax 0.8 \n")
draw.write("view ymin 0.85 \n")
draw.write("view ymax 1.25 \n")
draw.write("world xmin "+str(ax_min)+" \n")
draw.write("world xmax "+str(ax_max)+" \n")
draw.write("world ymin "+str(ax_min)+" \n")
draw.write("world ymax "+str(ay_max)+" \n")
draw.write("xaxis label \"Time [ns]\"  \n")
draw.write("xaxis label font \"Helvetica\"  \n")
draw.write("xaxis label char size 1.0 \n")
draw.write("xaxis ticklabel font \"Helvetica\"  \n")
draw.write("xaxis ticklabel char size 0.8 \n")
draw.write("xaxis tick major size 0.8 \n")
draw.write("xaxis tick minor size 0.4 \n")
draw.write("yaxis label \"RMSD [A]\"  \n")
draw.write("yaxis label font \"Helvetica\"  \n")
draw.write("yaxis label char size 1.0 \n")
draw.write("yaxis ticklabel font \"Helvetica\"  \n")
draw.write("yaxis ticklabel char size 0.8 \n")
draw.write("yaxis tick major size 0.8 \n")
draw.write("yaxis tick minor size 0.4 \n")
draw.write("legend 0.2, 0.7 \n")
draw.write("legend char size 0.8 \n")
draw.write("legend font \"Helvetica\" \n")
draw.close()

log = open("temp.log","w")
log.write(err)
log.close()