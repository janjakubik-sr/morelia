#This macro converts mol2 to yob runs energy minimization
OnError Exit

#scipt
DelAll
LoadMOL2 'target'
NameMol A,'MyName'
Style BallStick
Experiment Minimization
Experiment On
Wait ExpEnd
SaveYob 1,'output'
Clear

Exit
