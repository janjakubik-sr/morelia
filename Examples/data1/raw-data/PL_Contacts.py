# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import re
import numpy as np

#Input number of time-frames in trajectory
frames = 1000

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
                pattern = '_'+str(i)
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
plt.xticks(ind, residues)
plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0], p6[0]), ('HBond', 'Hydrophobic', 'Ion', 'Pi-Cation', 'PiPi', 'Water bridge'))

plt.show()
