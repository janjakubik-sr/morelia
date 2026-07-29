#This macro renames ligand for virtual screening
OnError Exit

#scipt
DelAll
LoadYob 'target'
NameMol A,'MyName'
SaveYob 1,'output'
Clear

Exit
