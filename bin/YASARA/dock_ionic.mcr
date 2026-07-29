# YASARA MACRO
# TOPIC:       Docking
# TITLE:       Analyzing ionic bonding of the ligand during a molecular dynamics simulation
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
  
# Create a table with 4 columns
MakeTab Final,Dimensions=2,Columns=4

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
  # Analyze ionic interactions
  resultlist() = ListIntMol (receptor),(ligand),Type=Ionic,Results=4
  if not count resultlist
    # If resultlist is empty, there are no ionic interactions
    Tabulate 0|
    Tabulate 0|
    Tabulate 'No ionic interactions found'
  else
    # Tabulate distance and the number of interactions
    Tabulate (resultlist3)|
    tabulate (resultlist4)|
    # Get a list of residues in receptor that are interacting with ligand via ionic interactions
    reslist() = ListIntRes (receptor),(ligand),Type=Ionic,Results=1
    # Convert list of unique residue IDs to residue names and numbers
    reslist() = ListRes (join reslist),Format=RES
    # Tabulate all residues as a single string
    Tabulate '(join reslist)'
  # Next pose
  i=i+1
  
# Save table, we got 4 columns
header='#____Pose| Distance|    Count| IBond Residues'
SaveTab Final,(MacroTarget)_ionic,Format=Text,Columns=4,NumFormat=9.2f,(header)

#HideMessage

# Exit YASARA if this macro was provided as command line argument in console mode
if runWithMacro and ConsoleMode
  Exit
