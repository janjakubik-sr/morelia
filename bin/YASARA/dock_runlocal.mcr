# YASARA MACRO
# TOPIC:       5. Structure prediction
# TITLE:       Docking a ligand to a receptor locally
# REQUIRES:    Structure
# AUTHOR:      Elmar Krieger
# LICENSE:     GPL
# DESCRIPTION: This macro samples local ligand conformations in a user-supplied complex. An analysis log file is written at the end.

# Parameter section - adjust as needed, but NOTE that some changes only take effect
# if you start an entirely new docking job, not if you continue an existing one. 
# =================================================================================

# You can either set the target structure by clicking on Options > Macro > Set target,
# or by uncommenting the line below and specifying it directly.
MacroTarget='macro_target'

# Docking method, either AutoDockLS or VINALS (LS=Local Search)
# Note that VINA local search yields only a single result. If you want an ensemble
# of ligand conformations, choose method 'AutoDockLS' and increase 'runs' below.
method = 'AutoDockLS'

# Number of docking runs. VINA local search always yields a single result only.
runs = number_of_runs

# Docking results usually cluster around certain hot spot conformations,
# and the lowest energy complex in each cluster is saved. Two complexes belong to
# different clusters if the ligand RMSD is larger than this minimum [A]:
rmsdmin = rmsd_min

# Set to 1 to keep the ligand completely rigid (alternatively you can provide
# the ligand as a *.yob file and fix certain dihedral angles only).
rigid = 0

# A selection of receptor residues to keep flexible, e.g. flexres = 'Lys 91 Leu 100'.
# Alternatively you can provide the receptor as a *.sce or *.yob file with fixed
# atoms, which gives better control (e.g. you can keep only part of a side-chain
# or even a terminal backbone flexible). 
flexres = '' 

# Normally no change required below this point
# ============================================
OnError Exit
# Sanity checks
if MacroTarget==''
  RaiseError "This macro requires a target. Either edit the macro file or click Options > Macro > Set target to choose a target structure"
if runs>999
  RaiseError "This macro can handle at most 999 docking runs"

structlist='receptor','ligand'

# Load receptor and ligand
Clear
# Do we already have a complex scene?
scefilename='(MacroTarget)_complex.sce'
scene = FileSize (scefilename)
if !scene
  # No, load PDB or YOB files of receptor and ligand.
  # Receptor and ligand coordinates are assumed to match
  for struct in structlist
    (struct)=0
    for type in 'yob','pdb'
      filename='(MacroTarget)_(struct).(type)'
      exists = FileSize (filename)
      if exists 
        if type=='pdb'
          (struct) = LoadPDB (filename),Center=No
        else
          (struct) = LoadYOb (filename)
        # Make sure coordinate systems match
        TransferObj ((struct)),(receptor)
        break
    if not (struct)
      RaiseError 'The (struct) was not found at (filename)'
  # Orient the scene, create a docking cell that includes the ligand. 
  NiceOriAll
  ZoomAll Steps=1
  Cell Auto,Extension=7,Obj (ligand)
  ForceField NOVA,SetPar=Yes
  Boundary Wall
  # Remember which atoms are fixed, and additionally fix the backbone
  fixedlist() = ListAtom fixed
  FixAtom Protein Backbone Obj (receptor)
  # Energy-minimize the receptor-ligand interaction
  Experiment Minimization
  Experiment On
  Wait ExpEnd
  Free
  if count fixedlist
    FixAtom (join fixedlist)
  # Save for reuse later
  SaveSce (scefilename)
else
  # Yes, a scene is present
  LoadSce (scefilename)
  # Verify that the cell is present
  simcell = CountObj SimCell
  if !simcell
    RaiseError 'If you provide a scene, it must contain a simulation cell, but none was found in (scefilename)'
  # The receptor is the object containing the first atom
  receptor = ListObj Atom 1
  # The ligand is the object containing the last atom
  ligand = ListObj Atom (Atoms)
  if ligand==receptor
    RaiseError 'Receptor and ligand must be present as separate objects'
# Force field used for charge assignment, not used by VINA.
ForceField AMBER03
# Docking is done without periodic boundaries
Boundary Wall
Longrange None
# Warn if the cell is too large
x,y,z = Cell
if x>96 or y>96 or z>96
  ShowMessage "A cell axis is longer than 96 A, which may reduce docking accuracy. Consider focusing on the active site."
  Wait ContinueButton
  HideMessage
# Rename receptor and ligand objects to make sure they can be identified later
NameObj (receptor),Receptor
NameObj (ligand),Ligand
# Keep selected side-chains flexible
if flexres!=''
  FreeObj Receptor
  FixAtom not Sidechain Res (flexres) and Obj Receptor 
# The segment name C*** is used to tag ligand Conformations, and cannot be used for the receptor
hit = ListAtom Segment C??? Obj Receptor
if hit
  MarkAtom (hit)
  segname = SegAtom (hit)
  RaiseError 'Segment names starting with "C" are unfortunately reserved to identify ligand conformations, please click "Edit > Rename > Segment" to rename the receptor segment "(segname)"' 
# The molecule name 'l' is reserved for the ligand
hit = ListAtom Mol l
if hit
  RaiseError 'The receptor must not contain a molecule named "l", please click Edit > Rename > Molecule and try again'
# Perform rigid docking if requested
if rigid
  FixAtom Obj Ligand
# Perform the docking with local sampling
Experiment Docking
  Method (method)
  ReceptorObj (receptor)
  LigandObj (ligand)
  Runs (runs)
  ClusterRMSD (rmsdmin)
  # Result file number must be separated with '_', in case MacroTarget also ends with number
  ResultFile (MacroTarget)_001
  # Uncomment below to set the number of energy evaluations (AutoDock ga_num_evals):
  # DockPar ga_num_evals 25000000
  # Uncomment the two lines below to provide your own atom parameters (VdW radii etc.)
  # GridPar parameter_file /Path/To/Custom/AD4_parameters.dat
  # DockPar parameter_file /Path/To/Custom/AD4_parameters.dat
  # If you want to keep all temporary files, uncomment below
  # TmpFileID YourChoice
Experiment On
Wait ExpEnd
# Save a scene with receptor and all ligand conformations
SaveSce (MacroTarget)

# Save a log file with an analysis
Console Off
RecordLog (MacroTarget)
print 'Local docking result analysis'
print '============================='
print
print '(runs) (method) docking runs of the ligand object (ligand) to the receptor object (receptor) yielded the following results,'
print 'sorted by binding energy: [More positive values indicates stronger binding, and negative values mean no binding]'
print
print 'Run |Bind.energy[kcal/mol]|Dissoc.constant[pM]| Contacting receptor residues'
print '----+---------------------+------------------+-----------------------------'
clusters=0
for i=001 to runs
  # Ligands have SegmentID C001, C002 etc.
  result = DockingResult i,'Obj (receptor)','Segment C(i)'
  print (result)
  # Count the clusters
  exists = FileSize (MacroTarget)_(i).yob
  if exists
    clusters=clusters+1
  if !(i%10)
    ShowMessage 'Analyzing docking results, (100*i/runs)% completed...'
    Wait 1
print
print 'After clustering the (runs) runs, the following (clusters) distinct complex conformations were found:'
print '[They all differ by at least (rmsdmin) A heavy atom RMSD]' 
print
print 'Cluster|Bind.energy[kcal/mol]|Dissoc. constant [pM]| Contacting receptor residues'
print '-------+---------------------+---------------------+-----------------------------'
for i=001 to clusters
  DelAll
  LoadYOb (MacroTarget)_(i)
  # Ligands have SegmentID C001, C002 etc.
  result = DockingResult i,'Segment !C???','Segment C???'
  print (result)
print 
print 'The results of the (runs) runs have been combined in a single scene at (MacroTarget).sce',Convert=No
print 'The (clusters) clusters have been saved as:'
for i=001 to clusters
  print '(MacroTarget)_(i).yob',Convert=No
StopLog
HideMessage
Console On
LoadSce (MacroTarget)
# Exit YASARA if this macro was provided as command line argument in console mode
if runWithMacro and ConsoleMode
  Exit

# EXTRACT DOCKING RESULT
# ======================
# Returns a result string. 'num' is a sequential number, 'recsel' selects the receptor, 'ligsel' the ligand.
def DockingResult num,recsel,ligsel
  # Get the binding energy from the B-factor...
  bindnrg = BFactorAtom (ligsel)
  # ...and the dissociation constant from the atomic property.
  dissconst = PropAtom (ligsel)
  if bindnrg<=0
    dissconst='               None'
  else
    dissconst= 00000000000000.0000+dissconst
  # Get receptor residues that contact the ligand
  reslist() = ListRes (recsel) with distance<4.0 from (ligsel),Format='RESNAME MOLNAME RESNUM'
  if !count reslist
    reslist='No residues in contact'
  else
    reslist=join reslist
  result='(000+num) |         (000000.0000+bindnrg) | (dissconst) | (reslist)'
  return result

