#This macro converts pdb to yasara version pdb 
OnError Exit

#Edit variables
MacroTarget='macro_target'

DelAll
LoadPDB (Macrotarget)
LoadPDB M2_3UON_OPM.pdb,Center=No,Correct=Yes
CleanAll
AlignObj 1,2,Method=MUSTANGPP
SavePDB 1,(Macrotarget),Format=PDB,Transform=Yes,UseCIF=No
Clear

Exit
