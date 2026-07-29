#This macro lists molecules in yob file
OnError Exit

#yob file
target = 'yob_file'

#scipt
DelAll
LoadYob (target).yob
ListMol Obj 1

Exit

