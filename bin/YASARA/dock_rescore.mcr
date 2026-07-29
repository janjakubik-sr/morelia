# YASARA MACRO
# TOPIC:       5. Structure prediction
# TITLE:       Rescore ligand poses with an energy minimization and local docking
# REQUIRES:    Structure
# AUTHOR:      Elmar Krieger
# LICENSE:     GPL
# DESCRIPTION: This macro loads a receptor-ligand complex obtained with YASARA's own docking macros or other software, energy-minimizes each ligand pose and rescores it using AutoDock VINA's local search, confined closely to the original ligand pose.

# Parameter section - adjust as needed, but NOTE that some changes only take effect
# if you start an entirely new docking job, not if you continue an existing one. 
# =================================================================================

# You can either set the target structure by clicking on Options > Macro > Set target,
# by providing it as command line argument (see docs at Essentials > The command line),
# or by uncommenting the line below and specifying it directly.
MacroTarget='macro_target'

# Docking method, either VINALS or AutoDockLS.
method='VINALS'

# Number of local docking runs per ligand (maximally 999, 0 lets YASARA choose based on the number of CPUs)
runs=0

# Sort results either by ligand number ('ligandnum'), binding energy ('bindnrg') or by
# efficiency ('effi'), i.e. binding energy per heavy atom, which addresses the bias towards
# larger ligands in screening.
sortby='bindnrg'

# Force field used to determine which ligand bonds can rotate, and also point charges if method='AutoDockLS'
forcefield='AMBER03'

# Selection of receptor atoms that are free to move during the minimization of each ligand pose.
# By default, only polar hydrogens are allowed to reorient.
recminsel='Element H with bond to Element N O P S'
# Alternatively, let all receptor side-chains move instead. Note that the more receptor atoms move, the better the
# ligand will bind. But since the energetic cost of the receptor's structural change is not considered in AutoDock's
# scoring function, moving receptor atoms may add noise instead of prediction accuracy.
#recminsel='Sidechain'
# Alternatively, let the entire receptor move.
# scoring function,  
#recminsel='All'
OnError Exit
Clear
# See if we are rescoring a previous normal YASARA docking run
scefilename='(MacroTarget).sce'
exists = FileSize (scefilename)
if exists
  # Use YASARA docking scene
  LoadSce (scefilename)
  receptor = ListObj Receptor
  ligandlist() = SplitObj Ligand
else
  # No, see if user provided a complex of receptor with ligands
  for type in 'yob','pdb'
    filename='(MacroTarget).(type)'
    exists = FileSize (filename)
    if exists 
      receptor = Load(type) (filename)
      break
  if !exists
    RaiseError 'Could neither find a YASARA scene (scefilename) from a previous docking, nor a complex in PDB or YOb format to import [e.g. (filename)]'
  # Check molecule names  
  recmolname = NameMol Atom 1
  ligmolname = NameMol Atom (Atoms)  
  if recmolname==ligmolname
    RaiseError 'The receptor and ligand poses share the same molecule name (recmolname), please click Edit > Name > Molecule and rename either of them'
  # Split the ligands into separate objects
  ligandlist() = SplitObj (receptor),Center=No,(ligmolname)
  DelVar ligandlist1
ligands = count ligandlist
# Check that ligands are present
if !ligands
  RaiseError 'No ligand poses found'
if !runs
  # For an optimum speed/accuracy tradeoff, we select at least 8 runs, but more if CPUs are available
  cpus = Processors
  runs=cpus
  while runs<8
    runs=runs+cpus
# Rename objects to make sure they can be identified later, ligands are 'Waiting' to be docked.
NameObj (receptor),Receptor
NameObj (join ligandlist),Waiting
RemoveObj (join ligandlist)
# Do we have a checkpoint file from a previous run?
checkpointfilename='(MacroTarget)_checkpoint.sce'
restarted = FileSize (checkpointfilename)
if restarted
  LoadSce (checkpointfilename)
Console Off
# Cycle through them
for i=1 to ligands
  ligand=ligandlist(i)
  ligname = NameObj (ligand)
  if ligname=='Done' or ligname=='Broken'
    continue
  AddObj (ligand)
  # Create a duplicate of the receptor to minimize (so that each minimization starts with the same receptor)
  recdupli = DuplicateObj (receptor)
  RemoveObj (receptor)
  # Energy-minimize the ligand and surrounding polar hydrogens
  FixObj (recdupli)
  FreeAtom (recminsel)
  Cell Auto,Extension=10
  ForceField NOVA,SetPar=Yes
  Boundary Periodic
  ShowMessage 'Energy-minimizing ligand (i)/(ligands)...'
  Experiment Minimization
  Experiment On
  Wait ExpEnd
  # Rotate the cell back to default orientation 
  alpha,beta,gamma = OriObj SimCell
  Rotate Y=(-beta)
  Rotate Z=(-gamma)
  Rotate X=(-alpha)
  # Use quick local docking to rescore the ligand
  ShowMessage 'Rescoring (i)/(ligands) with local docking...'
  # Create a cell tightly around the ligand by orienting the ligand along its major axes
  Cell Auto,Extension=2.5,Shape=Cuboid,Obj (ligand)
  Boundary Wall
  # Perform the docking with local sampling
  FreeObj (recdupli)
  ForceField (forcefield),SetPar=Yes
  Experiment Docking
    Method (method)
    ReceptorObj (recdupli)
    LigandObj (ligand)
    Runs (runs)
    # We only want a single result, the best complex
    ClusterRMSD 1000
    # Result file number must be separated with '_', in case MacroTarget also ends with number
    # If PDB files are needed too, replace .yob with .pdb to get both
    ResultFile (MacroTarget)_rescored_(000+i).yob 
  Experiment On
  Wait ExpEnd
  # Keep only the top scoring ligand pose
  DelMol Segment !C001 Obj (ligand)
  # Assign segment name for dock_play macro.
  segname = PoseSegName (i)
  SegObj (ligand),(segname)
  # Done, save checkpoint file
  NameObj (ligand),Done
  RemoveObj (ligand)
  # Bring back the original receptor for the next round
  DelObj (recdupli)
  AddObj (receptor)
  SaveSce (checkpointfilename)
# Join the ligands to a single object
AddObj (join ligandlist)
JoinObj (join ligandlist),(ligandlist1)  
NameObj (ligandlist1),Ligand
# Save the result scene for use by the dock_play macro
DelObj SimCell
SaveSce (MacroTarget)_rescored
# Save a log file with an analysis
RecordLog (MacroTarget)_rescored
AddObj all
print '#Ligand rescoring result'
print '#======================='
print '#'
print '#(ligands) ligand poses were rescored with an energy minimization followed by (runs) local docking runs each.'
print '#The following results were obtained sorted by binding energy:'
print '#[More positive energies indicate stronger binding, and negative energies mean no binding].'
print '#'
print "#Lig.|Effi[kcal/(mol*Atom)]|Bind.energy[kcal/mol]|Dissoc. constant [pM]| Name | Contacting receptor residues"
print '#----+---------------------+---------------------+---------------------+------+-----------------------------'
idx=1
for i=1 to ligands
  segname = PoseSegName (i)
  resultlist(idx),bindnrglist(idx),effilist(idx) = DockingResult i,'Obj (receptor)','Segment (segname) Obj (ligandlist1)'
  idx=idx+1
  if !(i%5)
    ShowMessage 'Analyzing docking results, (100*i/ligands)% completed...'
    Wait 1
if sortby!='ligandnum'
  # Sort results by binding energy or binding efficiency
  SortResults 'resultlist','(sortby)list'
for i=1 to count resultlist  
  print (resultlist(i))
print
rndseed = RandomSeed
print 'The random seed used during docking was (rndseed).'
print 'Point charges and dihedral barriers were obtained from the (forcefield) force field.'
StopLog
HideMessage
Console On

# Exit YASARA if this macro was provided as command line argument in console mode and not included from another macro
if runWithMacro and ConsoleMode and !IndentationLevel
  Exit

# EXTRACT DOCKING RESULT
# ======================
# Returns a result string and binding energy, extracted from atomic BFactor and Property.
# 'num' is a sequential number, 'recsel' selects the receptor, 'ligsel' the ligand.
def DockingResult num,recsel,ligsel
  # Get the binding energy from the B-factor...
  bindnrg = BFactorAtom (ligsel)
  # ...and the dissociation constant from the atomic property.
  dissconst = PropAtom (ligsel)
  if bindnrg<=0
    dissconst='               None'
  else
    dissconst= 00000000000000.0000+dissconst
  # Get compound name
  compound = CompoundMol (ligsel)
  if compound==''
    compound = NameRes (ligsel)
    compound=compound+' '
  # Get receptor residues that contact the ligand
  reslist() = ListRes (recsel) with distance<4.0 from (ligsel),Format='MOLNAME RESNAME RESNUM'
  if bindnrg<=-9999
    reslist='Broken molecule, no docking possible'
  elif !count reslist
    reslist='No residues in contact'
  else
    reslist=join reslist
  # Calculate ligand efficiency, i.e. binding energy divided by heavy atoms  
  heavyatoms = CountAtom Element !H (ligsel)
  effi = bindnrg/heavyatoms
  result='(0000+num)|         (000000.0000+effi) |         (000000.0000+bindnrg) | (dissconst) | (compound) | (reslist)'
  return result,bindnrg,effi

# SORT KEYLIST BY VALUES, HIGHEST FIRST
# =====================================
# Lists are passed by reference, see Yanaconda doc section 'Calls to user-defined functions..'
def SortResults keylist,vallist
  caller (keylist),(vallist)

  elements=count (keylist)
  do
    sorted=1
    for i=1 to elements-1
      if (vallist)(i)<(vallist)(i+1)
        swap=(vallist)(i)
        (vallist)(i)=(vallist)(i+1)
        (vallist)(i+1)=swap
        swap=(keylist)(i)
        (keylist)(i)=(keylist)(i+1)
        (keylist)(i+1)=swap
        sorted=0
  while not sorted

# DETERMINE LIGAND POSE SEGMENT NAME
# ==================================
def PoseSegName number
  if number<1000
    name='C(000+number)'
  else
    name='(0000+number)'
  return name
