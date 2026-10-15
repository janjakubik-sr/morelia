#This macro converts smi to yob runs energy minimization
OnError Exit

#scipt
DelAll
LoadSMI 'target'
NameMol A,'MyName'
Style BallStick
Experiment Minimization
Experiment On
Wait ExpEnd
SaveYob 1,'output'
Clear

Exit
