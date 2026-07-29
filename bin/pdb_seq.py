#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 14:15:28 2023

@author: Jan Jakubik jan.jakubik@fgu.cas.cz
"""
import os
import time
import wx

from Bio.PDB import PDBParser
from Bio.PDB import PPBuilder
parser = PDBParser(PERMISSIVE=1)
ppb = PPBuilder()

#error log
if os.path.isfile("temp.log"):
    os.remove("temp.log")
err = ''

#select pdf file
dialog = wx.FileDialog(None, "Choose pdb file", os.getcwd(), "","*.pdb", wx.FD_OPEN)
if dialog.ShowModal() == wx.ID_OK:
    pdb=dialog.GetPath()
    dialog.Destroy()
    #get base name
    full_path=os.path.splitext(pdb)[0]
    structure_id=os.path.basename(full_path).split('/')[-1]
    structure = parser.get_structure(structure_id, pdb)
    #prepare file for writing results
    if os.path.isfile("pdb_seq.res"):
        os.remove("pdb_seq.res")
    res = open("pdb_seq.res","w")
    res.write(str(structure_id)+".pdb description\n")
    for model in structure:
        res.write('Model '+str(model.id)+':\n')
        for chain in model:
            if len(chain) == 1:
                res.write('Chain '+str(chain.id)+' is hetero residue\n\n')
            else:
                res.write('Chain ' +str(chain.id)+': '+str(len(chain))+' residues\n')
                for pp in ppb.build_peptides(chain):
                    res.write(str(pp.get_sequence())+'\n\n')
    now = time.strftime("%d.%m.%Y %H:%M:%S")
    res.write("\n\n"+str(now)+"\n\n")
    res.close()
else:
    msg = ('No file selected!\nExiting!')
    err = msg
    dialog = wx.MessageDialog(None, msg, 'Error', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()

log = open("temp.log","w")
log.write(err)
log.close()