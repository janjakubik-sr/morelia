#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 15:06:35 2023

@author: roshi
"""

from Bio.PDB import PDBParser
from Bio.PDB import PPBuilder
from Bio.PDB import DSSP
parser = PDBParser(PERMISSIVE=0)
ppb = PPBuilder()

structure_id = "3uon"
filename = "3uon.pdb"
structure = parser.get_structure(structure_id, filename)

for model in structure:
    print('Model '+str(model.id)+':')
    for chain in model:
        if len(chain) == 1:
            print('Chain '+str(chain.id)+' is hetero residue')
        else:
            print('Chain ' +str(chain.id)+': '+str(len(chain))+' residues')
            for pp in ppb.build_peptides(chain):
                dssp = DSSP(chain, filename)
                print(pp.get_sequence())