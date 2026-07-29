#This macro lists HBonds between ligand and receptor in yob file
OnError Exit

#yob file
target = 'yob_file'

# Receptor and ligand molecules
receptor = 'mol_receptor'
ligand = 'mol_ligand'

#scipt
DelAll
LoadYob (target).yob

# Create a table with 5 columns
MakeTab Final,Dimensions=2,Columns=5

resultlist() = ListHBoRes (receptor),(ligand),Results=4
if not count resultlist
  # If resultlist is empty, there are no hydrogen bonds
  Tabulate 0|
  Tabulate 0|
  Tabulate 0|
  Tabulate 'No hydrogen bonds found'
else
  # Tabulate energy and the numbers of accepted and donated bonds
  Tabulate (resultlist2)|
  Tabulate (resultlist3)|
  Tabulate (resultlist4)|
  # Get a list of residues in receptor that are H-bonded to ligand
  reslist() = ListHBoRes (receptor),(ligand),Results=1
  # Convert list of unique residue IDs to residue names and numbers
  reslist() = ListRes (join reslist),Format=RES
  # Tabulate all residues as a single string
  Tabulate '(join reslist)'

# Save table, we got 5 columns
header='#____Pose|   Energy| Accepted|  Donated| HBond Residues'
SaveTab Final,HBonds,Format=Text,Columns=5,NumFormat=9.2f,(header)

Exit

