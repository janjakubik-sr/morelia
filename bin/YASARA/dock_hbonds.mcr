# YASARA MACRO
# TOPIC:       Docking
# TITLE:       Analyzing hydrogen bonding of the ligand from docking
# REQUIRES:    Dynamics
# AUTHOR:      Jan Jakubik
# LICENSE:     GPL
# DESCRIPTION: This macro analyzes hydrogen bonds between ligand and receptor of individual docking poses

# You can either set the target structure by clicking on Options > Macro > Set target,
# or by uncommenting the line below and specifying it directly.
MacroTarget = 'macro_target'

# Forcefield to use (these are all YASARA commands, so no '=' used)
ForceField AMBER03,SetPar=Yes

# Receptor and ligand molecules
receptor = 'mol_receptor'
ligand = 'mol_ligand'

# First and last pose to be analysed
firstpose=first_pose
lastpose=last_pose

# No change required below this point
# ===================================
OnError Exit
# Do we have a target?
if MacroTarget==''
  RaiseError "This macro requires a target. Either edit the macro file or click Options > Macro > Set target to choose a target structure"
else
  LoadSce (MacroTarget)
  
# Create a table with 5 columns
MakeTab Final,Dimensions=2,Columns=5

#Clear
#Console Off

i=000+firstpose
last=000+lastpose
while i<=last
  DelAll
  # Load next pose
  LoadYob (MacroTarget)_(i)
  # Add pose number to the table.
  ShowMessage 'Analyzing pose (0+i)'
  Tabulate (0+i)|
  Wait 1
  # Analyze hydrogen bonds
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
  # Next pose
  i=i+1
  
# Save table, we got 5 columns
header='#____Pose|   Energy| Accepted|  Donated| HBond Residues'
SaveTab Final,(MacroTarget)_hbonds,Format=Text,Columns=5,NumFormat=9.2f,(header)

#HideMessage

# Exit YASARA if this macro was provided as command line argument in console mode
if runWithMacro and ConsoleMode
  Exit
