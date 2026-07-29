#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import re
import wx
import subprocess
import numpy as np
import preferences as pref
install_dir = (pref.install_dir)

#error log
if os.path.isfile('temp.log'):
    os.remove('temp.log')
err = ''

#Input number of time-frames in trajectory
ask = ('Enter number of frames\n')
dlg = wx.TextEntryDialog(None, ask,"Input","1000")
if dlg.ShowModal() == wx.ID_OK:
    frames = dlg.GetValue()
    dlg.Destroy()
    frames = float(frames)
else:
    msg = ('Aborted!\nExiting!')
    err = msg
    log = open("temp.log","w")
    log.write(err)
    log.close()
    dialog = wx.MessageDialog(None, msg, 'Error', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()
    exit()

#Extract data from individual files
HB_res = []
with open('PL-Contacts_HBond.dat') as HB:
    for line in HB:
        fields = line.strip().split()
        HB_res.append(fields[1])

HP_res = []
with open('PL-Contacts_Hydrophobic.dat') as HP:
    for line in HP:
        fields = line.strip().split()
        HP_res.append(fields[1])

Ion_res = []
with open('PL-Contacts_Ionic.dat') as Ion:
    for line in Ion:
        fields = line.strip().split()
        Ion_res.append(fields[1])

PiCat_res = []
with open('PL-Contacts_Pi-Cation.dat') as PiCat:
    for line in PiCat:
        fields = line.strip().split()
        PiCat_res.append(fields[1])

PiPi_res = []
with open('PL-Contacts_Pi-Pi.dat') as PiPi:
    for line in PiPi:
        fields = line.strip().split()
        PiPi_res.append(fields[1])

WB_res = []
with open('PL-Contacts_WaterBridge.dat') as WB:
    for line in WB:
        fields = line.strip().split()
        WB_res.append(fields[1])

#Count numer of intractions of individual types
sum_data = np.array([0,0,0,0,0,0,0])
for i in range(1, 999):
    n_HB_res = HB_res.count(str(i))
    n_HP_res = HP_res.count(str(i))
    n_Ion_res = Ion_res.count(str(i))
    n_PiCat_res = PiCat_res.count(str(i))
    n_PiPi_res = PiPi_res.count(str(i))
    n_WB_res = WB_res.count(str(i))
    total = n_HB_res + n_HP_res + n_Ion_res + n_PiCat_res + n_PiPi_res + n_WB_res
    total_min = 0.05 * frames
    if total > total_min:
        with open('P_RMSF.dat') as ID:
            for line in ID:
                pattern = '_'+str(i)+' '
                if re.search(pattern, line):
                    res_id = (line.strip().split())[2]
        exec('sum_data_'+str(i)+'= np.stack((res_id, n_HB_res, n_HP_res, n_Ion_res, n_PiCat_res, n_PiPi_res, n_WB_res), axis=-1)')
        exec('sum_data = np.vstack((sum_data,sum_data_'+str(i)+'))')

#Save summary data to file
sum_data = np.delete(sum_data,0,0)
with open('PL_Contacts_Summary.dat','w') as res:
    np.savetxt(res,sum_data,fmt='%s')

#Compute fractions of interactions
residues = []
HB_count = []
HP_count = []
Ion_count = []
PiCat_count = []
PiPi_count = []
WB_count = []
with open('PL_Contacts_Summary.dat') as SD:
    for line in SD:
        fields = line.strip().split()
        residues.append(fields[0])
        HB_count.append(float(fields[1])/frames)
        HP_count.append(float(fields[2])/frames)
        Ion_count.append(float(fields[3])/frames)
        PiCat_count.append(float(fields[4])/frames)
        PiPi_count.append(float(fields[5])/frames)
        WB_count.append(float(fields[6])/frames)


#Plot with matplotlib
import matplotlib.pyplot as plt

N = len(residues)
ind = np.arange(N)
b3 = np.add(HB_count,HP_count)
b4 = np.add(b3,Ion_count)
b5 = np.add(b4,PiCat_count)
b6 = np.add(b5,PiPi_count)
p1 = plt.bar(ind, HB_count)
p2 = plt.bar(ind, HP_count, bottom = HB_count)
p3 = plt.bar(ind, Ion_count, bottom = b3)
p4 = plt.bar(ind, PiCat_count, bottom = b4)
p5 = plt.bar(ind, PiPi_count, bottom = b5)
p6 = plt.bar(ind, WB_count, bottom = b6)

plt.xlabel('Residue #')
plt.ylabel('Interactiions fraction')
plt.xticks(ind, residues, rotation=45)
plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0], p6[0]), ('HBond', 'Hydrophobic', 'Ion', 'Pi-Cation', 'PiPi', 'Water bridge'))

plt.show()

#Plot with Grace

s0 = np.column_stack((ind, HB_count))
np.savetxt('s0.dat',s0,fmt='%s')
s1 = np.column_stack((ind, HP_count))
np.savetxt('s1.dat',s1,fmt='%s')
s2 = np.column_stack((ind, Ion_count))
np.savetxt('s2.dat',s2,fmt='%s')
s3 = np.column_stack((ind, PiCat_count))
np.savetxt('s3.dat',s3,fmt='%s')
s4 = np.column_stack((ind, PiPi_count))
np.savetxt('s4.dat',s4,fmt='%s')
s5 = np.column_stack((ind, WB_count))
np.savetxt('s5.dat',s5,fmt='%s')

#prepare file for Grace plot
if os.path.isfile("PL_Contacts.draw"):
    os.remove("PL_Contacts.draw")
draw = open("PL_Contacts.draw","w")

ax_min = -1
ax_max = N
ay_min = 0
ay_max = 2

draw.write("page layout portrait \n")
draw.write("view xmin 0.2 \n")
draw.write("view xmax 0.8 \n")
draw.write("view ymin 0.85 \n")
draw.write("view ymax 1.25 \n")

draw.write(" G0  type Chart \n")
draw.write(" G0  stacked true \n")
draw.write("  xaxis  label \"Residue \" \n")
draw.write("  xaxis  label font \"Helvetica\" \n")
draw.write("  xaxis  label char size 1.0 \n")
draw.write("  xaxis  ticklabel font \"Helvetica\" \n")
draw.write("  xaxis  ticklabel char size 0.8 \n")
draw.write("  xaxis  ticklabel angle 45 \n")
draw.write("  yaxis  label \"Interactions fraction\" \n")
draw.write("  yaxis  label font \"Helvetica\" \n")
draw.write("  yaxis  label char size 1.0 \n")
draw.write("  yaxis  ticklabel font \"Helvetica\" \n")
draw.write("  yaxis  ticklabel char size 0.8 \n")
draw.write("  Legend  anchor 0.22 , 0.64 \n")
draw.write("  Legend  font \"Helvetica\" \n")
draw.write("  xaxis  tick off \n")
draw.write("  xaxis  tick spec type both \n")
draw.write("  xaxis  tick spec "+str(N)+" \n")
for j in ind:
    draw.write("xaxis  tick major "+str(j)+" ,"+str(j)+" \n")
    draw.write("xaxis  ticklabel "+str(j)+" ,\""+str(residues[j])+"\" \n")

draw.write(" read xy \"s0.dat\" \n")
draw.write("    S0 type bar \n")
draw.write("    S0 symbol size 1.500000 \n")
draw.write("    S0 symbol color 22 \n")
draw.write("    S0 symbol fill color 22 \n")
draw.write("    S0 symbol fill pattern 1 \n")
draw.write("    S0 symbol linestyle 0 \n")
draw.write("    S0 line type 0 \n")
draw.write("    S0 legend on \n")
draw.write("    S0 legend \"HBond\" \n")

draw.write(" read xy \"s1.dat\" \n")
draw.write("    S1 type bar \n")
draw.write("    S1 symbol size 1.500000 \n")
draw.write("    S1 symbol color 28 \n")
draw.write("    S1 symbol fill color 28 \n")
draw.write("    S1 symbol fill pattern 1 \n")
draw.write("    S1 symbol linestyle 0 \n")
draw.write("    S1 line type 0 \n")
draw.write("    S1 legend on \n")
draw.write("    S1 legend \"Hydrophobic\" \n")

draw.write(" read xy \"s2.dat\" \n")
draw.write("    S2 type bar \n")
draw.write("    S2 symbol size 1.500000 \n")
draw.write("    S2 symbol color 32 \n")
draw.write("    S2 symbol fill color 32 \n")
draw.write("    S2 symbol fill pattern 1 \n")
draw.write("    S2 symbol linestyle 0 \n")
draw.write("    S2 line type 0 \n")
draw.write("    S2 legend on \n")
draw.write("    S2 legend \"Ionic\" \n")

draw.write(" read xy \"s3.dat\" \n")
draw.write("    S3 type bar \n")
draw.write("    S3 symbol size 1.500000 \n")
draw.write("    S3 symbol color 24 \n")
draw.write("    S3 symbol fill color 24 \n")
draw.write("    S3 symbol fill pattern 1 \n")
draw.write("    S3 symbol linestyle 0 \n")
draw.write("    S3 line type 0 \n")
draw.write("    S3 legend on \n")
draw.write("    S3 legend \"Pi-Cation\" \n")

draw.write(" read xy \"s4.dat\" \n")
draw.write("    S4 type bar \n")
draw.write("    S4 symbol size 1.500000 \n")
draw.write("    S4 symbol color 16 \n")
draw.write("    S4 symbol fill color 16 \n")
draw.write("    S4 symbol fill pattern 1 \n")
draw.write("    S4 symbol linestyle 0 \n")
draw.write("    S4 line type 0 \n")
draw.write("    S4 legend on \n")
draw.write("    S4 legend \"Pi-Pi\" \n")

draw.write(" read xy \"s5.dat\" \n")
draw.write("    S5 type bar \n")
draw.write("    S5 symbol size 1.500000 \n")
draw.write("    S5 symbol color 41 \n")
draw.write("    S5 symbol fill color 41 \n")
draw.write("    S5 symbol fill pattern 1 \n")
draw.write("    S5 symbol linestyle 0 \n")
draw.write("    S5 line type 0 \n")
draw.write("    S5 legend on \n")
draw.write("    S5 legend \"Water bridge\" \n")

draw.write("world ymin "+str(ay_min)+" \n")
draw.write("world ymax "+str(ay_max)+" \n")
draw.write("world xmin "+str(ax_min)+" \n")
draw.write("world xmax "+str(ax_max)+" \n")
draw.write("  yaxis  tick major 0.5 \n")
draw.write("  yaxis  tick minor ticks 1 \n")

draw.close()

cmd = (str(install_dir)+'3_to_1_letter.sh')
p1 = subprocess.Popen(cmd, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
p1.wait()

log = open("temp.log","w")
log.write(err)
log.close()
