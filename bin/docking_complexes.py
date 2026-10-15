# -*- coding: utf-8 -*-
"""
docking_complexes.py -- create YASARA docking complex scenes.

In the working directory (e.g. dock/M1) the script expects a *_complex.sce
file (e.g. M1_5CXV_ortho_complex.sce) that contains the receptor and the
docking box (simulation cell), and ligands as .yob files. For each ligand,
the ligand is added to the complex scene and the scene is saved as
receptor-base-name_ligand-name_complex.sce into a subfolder named
receptor-base-name_ligand-name (e.g. M1_5CXV_ortho_JRPM-01).
"""

#No edits below

import os
import glob
import wx
import subprocess
import preferences as pref
install_dir = (pref.install_dir)

#error log
if os.path.isfile("temp.log"):
    os.remove("temp.log")
err = ''

#select sce file
dialog = wx.FileDialog(None, "Choose receptor complex sce file", os.getcwd(), "", "*_complex.sce", wx.FD_OPEN)
if dialog.ShowModal() == wx.ID_OK:
    sce = dialog.GetPath()
    dialog.Destroy()
    workdir = os.path.dirname(os.path.abspath(sce))
    os.chdir(workdir)
    #get receptor base name (strip trailing _complex)
    target = os.path.splitext(os.path.basename(sce))[0].replace("_complex", "")
    #find ligands (yob files in working directory)
    ligands = sorted(os.path.splitext(os.path.basename(f))[0]
                     for f in glob.glob("*.yob"))
    if not ligands:
        err = "No yob ligands found."
        dialog = wx.MessageDialog(None, 'No ligand .yob files found in\n'+workdir+'\nExiting!', 'Error', wx.OK)
        dialog.ShowModal()
        dialog.Destroy()
    else:
        msg = ('Complex scenes will be created for '+str(len(ligands))+' ligands:\n\n'+"\n".join(ligands))
        dialog = wx.MessageDialog(None, msg, 'Info', wx.OK | wx.CANCEL)
        proceed = dialog.ShowModal()
        dialog.Destroy()
        if proceed == wx.ID_OK:
            #create per-ligand subfolders
            for ligand in ligands:
                outdir = os.path.join(target+"_"+ligand)
                if not os.path.isdir(outdir):
                    os.makedirs(outdir)
            #prepare macro
            liglist = ",".join("'"+lig+"'" for lig in ligands)
            with open(os.path.join(str(install_dir), 'YASARA', 'docking_complex.mcr')) as f:
                macro = f.read()
            macro = macro.replace('macro_target', str(target))
            macro = macro.replace('ligand_list', liglist)
            with open('docking_complex.mcr', 'w') as f:
                f.write(macro)
            #run macro
            cmd3 = ('yasara -txt docking_complex.mcr > docking_complexes.log')
            p3 = subprocess.Popen(cmd3, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
            dialog = wx.MessageDialog(None, 'Creating docking complexes\nPlease wait', 'Please wait', wx.OK)
            dialog.ShowModal()
            p3.wait()
            dialog.Destroy()
            #check results
            done = [lig for lig in ligands
                    if os.path.isfile(os.path.join(target+"_"+lig, target+"_"+lig+"_complex.sce"))]
            msg = ('Created complex scenes for '+str(len(done))+' of '+str(len(ligands))+' ligands.')
            dialog = wx.MessageDialog(None, msg, 'Done', wx.OK)
            dialog.ShowModal()
            dialog.Destroy()
            err = msg
        else:
            err = "Canceled"
else:
    msg = ('No sce input file selected!\nExiting!')
    err = "Canceled"
    dialog = wx.MessageDialog(None, msg, 'Error', wx.OK)
    dialog.ShowModal()
    dialog.Destroy()

log = open("temp.log", "w")
log.write(err)
log.close()
