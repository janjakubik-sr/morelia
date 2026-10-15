# YASARA MACRO
# TOPIC:       5. Structure prediction
# TITLE:       Create receptor-ligand complex scenes for docking
# REQUIRES:    Structure
# DESCRIPTION: Loads a receptor complex scene (receptor + docking box), adds each
#              ligand YOb to it and saves one complex scene per ligand into a
#              subfolder named (MacroTarget)_(ligand).

# Parameter section - adjusted by docking_complexes.py
# =====================================================

# Base name of the receptor complex scene, '(MacroTarget)_complex.sce'
MacroTarget='macro_target'

# List of ligand YOb files (without extension)
liglist = ligand_list

# Normally no change required below this point
# ============================================
OnError Exit
Clear
# Load the receptor complex scene and verify that the cell is present
scefilename='(MacroTarget)_complex.sce'
exists = FileSize (scefilename)
if !exists
  RaiseError 'Complex scene (scefilename) not found'
for ligand in liglist
  LoadSce (scefilename)
  simcell = CountObj SimCell
  if !simcell
    RaiseError 'The scene (scefilename) must contain a simulation cell (docking box), but none was found'
  # The receptor is the object containing the first atom
  receptor = ListObj Atom 1
  lig = LoadYOb (ligand)
  if !lig
    RaiseError 'The ligand was not found at (ligand).yob'
  # Make sure coordinate systems match
  TransferObj (lig),(receptor)
  # Rename objects so that docking macros can identify them later
  NameObj (receptor),Receptor
  NameObj (lig),Ligand
  SaveSce '(MacroTarget)_(ligand)/(MacroTarget)_(ligand)_complex.sce'
  DelObj (lig)
# Exit YASARA if this macro was provided as command line argument in console mode
if runWithMacro and ConsoleMode
  Exit
