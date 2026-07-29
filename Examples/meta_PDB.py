#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 17:04:44 2023

@author: roshi
"""

from Bio.PDB import PDBParser
from Bio.PDB import PPBuilder
parser = PDBParser(PERMISSIVE=1)
ppb = PPBuilder()

structure_id = "6oik"
filename = "6oik.pdb"
structure = parser.get_structure(structure_id, filename)
print('Structure '+str(structure.id)+':')
for model in structure:
    print('Model '+str(model.id)+':')
    for chain in model:
        print('Chain ' +str(chain.id))
        if len(chain) == 1:
            print(' has 1 residue')
        else:
            print(' has '+str(len(chain))+' residues')
