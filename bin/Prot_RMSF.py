# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#No edits below

import os
import sys
import wx
import numpy as np

#Extract data
frame = []
contact = []
rmsf = []
with open('P_RMSF.dat') as RMSF:
    for line in RMSF:
        fields = line.strip().split()
        frame.append(fields[0])
        rmsf.append(fields[5])
        if fields[3] == 'Yes':
            contact.append(fields[5])
        else:
            contact.append('0')

Prot_RMSF = np.column_stack((frame, rmsf))
Prot_RMSF = np.delete(Prot_RMSF,0,0)
np.savetxt('Prot_RMSF.dat',Prot_RMSF,fmt='%s')

PL_contacts = np.column_stack((frame, contact))
PL_contacts = np.delete(PL_contacts,0,0)
np.savetxt('PL_contacts.dat',PL_contacts,fmt='%s')

data = np.loadtxt('Prot_RMSF.dat')

#Plot with matplotlib
import matplotlib.pyplot as plt
x_max = np.amax(data[:,0]) * 1.1
y_max = np.amax(data[:,1]) * 1.1
plt.title("Protein RMSF")
plt.axis([0,x_max,0,y_max])
plt.plot(data[:,0],data[:,1])
plt.xlabel('Residue #')
plt.ylabel('RMSF [A]')
plt.show()

#Plot with Grace
#prepare file for Grace plot
if os.path.isfile("Prot_RMSF.draw"):
    os.remove("Prot_RMSF.draw")
draw = open("Prot_RMSF.draw","w")

ax_min = 0
ax_max = x_max
ay_min = 0
ay_max = y_max

draw.write("read xy \"Prot_RMSF.dat\" \n")
draw.write("s0 symbol 0 \n")
draw.write("s0 line type 1 \n")
draw.write("s0 line color 1 \n")

draw.write("read xy \"PL_contacts.dat\" \n")
draw.write("s1 symbol 0 \n")
draw.write("s1 line type 0 \n")
draw.write("s1 line color 2 \n")
draw.write("s1 dropline on \n")

draw.write("page layout portrait \n")
draw.write("view xmin 0.2 \n")
draw.write("view xmax 0.8 \n")
draw.write("view ymin 0.85 \n")
draw.write("view ymax 1.25 \n")
draw.write("world xmin "+str(ax_min)+" \n")
draw.write("world xmax "+str(ax_max)+" \n")
draw.write("world ymin "+str(ax_min)+" \n")
draw.write("world ymax "+str(ay_max)+" \n")
draw.write("xaxis label \"Residue #\"  \n")
draw.write("xaxis label font \"Helvetica\"  \n")
draw.write("xaxis label char size 1.0 \n")
draw.write("xaxis ticklabel font \"Helvetica\"  \n")
draw.write("xaxis ticklabel char size 0.8 \n")
draw.write("xaxis tick major size 0.8 \n")
draw.write("xaxis tick minor size 0.4 \n")
draw.write("yaxis label \"RMSF [A]\"  \n")
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
