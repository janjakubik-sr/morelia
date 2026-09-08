#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 16:01:00 2016

@author: Jan Jakubik jan.jakubik@fgu.cas.cz
"""
#!/usrenv python
import wx
import os
import sys
import subprocess
import time

import preferences as pref
install_dir = (pref.install_dir)
sys.path.append(install_dir)
doc_dir = (pref.doc_dir)
sys.path.append(doc_dir)
example_dir = (pref.example_dir)
sys.path.append(example_dir)
width = (pref.width)
height = (pref.height)

#temp file
if os.path.isfile("nc.tmp"):
    os.remove("nc.tmp")

def exec_full(filepath):
    global_namespace = {
        "__file__": filepath,
        "__name__": "__main__",
    }
    with open(filepath, 'rb') as file:
        exec(compile(file.read(), filepath, 'exec'), global_namespace)


welcome = "\n******************************************\nHello, welcome to Morelia\nMolecular modelling utilities written in python\nby Jan Jakubik jan.jakubik@fgu.cas.cz\nCreative Commons Licence BY-NC\n*****************************************\n"

stamp = time.strftime("%Y%m%d%H%M%S")
class MainWindow(wx.Frame):
    def __init__(self, filename=(stamp+'.log')):
        super(MainWindow, self).__init__(None, size=(width,height))
        self.filename = filename
        self.dirname = '.'
        self.CreateInteriorWindowComponents()
        self.CreateExteriorWindowComponents()

    def CreateInteriorWindowComponents(self):
        self.logger = wx.TextCtrl(self, value=welcome, style=wx.TE_MULTILINE|wx.TE_READONLY)
        
    def CreateExteriorWindowComponents(self):
        self.CreateMenu()
        self.CreateStatusBar()
        self.SetTitle()

    def CreateMenu(self):
        fileMenu = wx.Menu()
        for id, label, helpText, handler in \
            [(wx.ID_SAVE, '&Save log', 'Save the current log file', self.OnSave),
             (wx.ID_OPEN, '&Working Directory', 'Change working directory', self.OnCWD),
             (wx.ID_PREFERENCES, 'Preferences', 'Define default aplications, etc.', self.OnPrefs),
             (wx.ID_EXIT, 'E&xit', 'Terminate the program', self.OnExit)
            ]:
            if id == None:
                fileMenu.AppendSeparator()
            else:
                item = fileMenu.Append(id, label, helpText)
                self.Bind(wx.EVT_MENU, handler, item)

        autoDockMenu = wx.Menu()
        for id, label, helpText, handler in \
            [(wx.ID_ANY, 'smi to pdbqt', 'Convert SMILES files to PDBQT', self.OnSmi2Pdbqt),
             (wx.ID_ANY, 'pdb to pdbqt', 'Prepare receptor PDB to PDBQT', self.OnPdb2Pdbqt),
             (wx.ID_ANY, 'Map receptor', 'Run AutoGrid on receptor GPF', self.OnMapReceptor),
             (wx.ID_ANY, 'Batch docking', 'Dock ligands to receptor with AutoDock-GPU', self.OnBatchDock),
	      ]:
            if id == None:
                autoDockMenu.AppendSeparator()
            else:
                item = autoDockMenu.Append(id, label, helpText)
                self.Bind(wx.EVT_MENU, handler, item)

        DesmondMenu = wx.Menu()
        for id, label, helpText, handler in \
            [(wx.ID_ANY, 'Align Trajectory', 'Align Trajectory', self.OnTrjAlign),
             (wx.ID_ANY, 'Ligand RMSD', 'Ligand RMSD', self.OnLigRMSD),
             (wx.ID_ANY, 'Ligand Contact Histogram', 'Protein-Ligand Contact Histogram', self.OnPLCon),
             (wx.ID_ANY, 'Minimization', 'Energy Minimization', self.OnMMin),
             (wx.ID_ANY, 'Molecular Dynamics', 'Molecular Dynamics', self.OnMD),
             (wx.ID_ANY, 'MD dummy host', 'MD dummy host', self.OnDummy),
             (wx.ID_ANY, 'Protein backbone RMSD', 'Protein backbone RMSD', self.OnPBBRMSD),
             (wx.ID_ANY, 'Protein RMSF', 'Protein RMSF', self.OnPRMSF),
             (wx.ID_ANY, 'Simulation Interaction Diagram', 'Schrodinger Simulation Interaction Diagram', self.OnSID),
             (wx.ID_ANY, 'Simulation Quality Analysis', 'Simulation Quality Analysis', self.OnQA),
             (wx.ID_ANY, 'Trajectory RMSD Clustering', 'Trajectory RMSD Clustering', self.OnTrjRMSD),
	      ]:
            if id == None:
                DesmondMenu.AppendSeparator()
            else:
                item = DesmondMenu.Append(id, label, helpText)
                self.Bind(wx.EVT_MENU, handler, item)

        GromacsMenu = wx.Menu()
        for id, label, helpText, handler in \
            [(wx.ID_ANY, 'Compress .trr to .xtc', 'Compress GROMACS .trr trajectory to .xtc', self.OnTRR2XTC),
             (wx.ID_ANY, 'Desmond MM-PBSA', 'Run Desmond gmx_MMPBSA workflow', self.OnDesmondMMPBSA),
	      ]:
            if id == None:
                GromacsMenu.AppendSeparator()
            else:
                item = GromacsMenu.Append(id, label, helpText)
                self.Bind(wx.EVT_MENU, handler, item)

        pyMenu = wx.Menu()
        for id, label, helpText, handler in \
            [(wx.ID_ANY, 'Measure ligand atoms distance to receptor AA centroid', 'Measure distance', self.OnLigDistance),
	      ]:
            if id == None:
                pyMenu.AppendSeparator()
            else:
                item = pyMenu.Append(id, label, helpText)
                self.Bind(wx.EVT_MENU, handler, item)

        vmdMenu = wx.Menu()
        for id, label, helpText, handler in \
            [(wx.ID_ANY, 'Convert .cms to .gro', 'Convert Desmond .cms structure to GROMACS .gro', self.OnCMS2GRO),
             (wx.ID_ANY, 'Convert .dtr to .trr', 'Convert Desmond .dtr trajectory to GROMACS .trr', self.OnDTR2TRR),
	      ]:
            if id == None:
                vmdMenu.AppendSeparator()
            else:
                item = vmdMenu.Append(id, label, helpText)
                self.Bind(wx.EVT_MENU, handler, item)

        yasaraMenu = wx.Menu()
        for id, label, helpText, handler in \
            [(wx.ID_ANY, 'Align PDB to membrane', 'PDB to OPM', self.OnPDBOPM),
             (wx.ID_ANY, 'Archive current folder', 'Compress to archive and clean current folder', self.OnArchive),
             (wx.ID_ANY, 'Dock - Local search', 'Dock using AutoLS', self.OnDockLS),
             (wx.ID_ANY, 'Dock - Rescore', 'Rescore docking poses', self.OnDockRe),
             (wx.ID_ANY, 'Dock - Extract top poses', 'Extract top poses from docking archives', self.OnExtractTop),
             (wx.ID_ANY, 'Liagnd-Receptor interactions', 'Analyze interactions', self.OnInteractions),
             (wx.ID_ANY, 'Make Si ligands topologies', 'Generate YASARA usertopo for Si ligands', self.OnMakeSiTopo),
             (wx.ID_ANY, 'Minimize Si ligands', 'Minimize Si ligands with YASARA', self.OnMinSiLigands),
             (wx.ID_ANY, 'Hydrogen bonds', 'Find hydrogen bonds in the docking poses', self.OnHBonds),
             (wx.ID_ANY, 'Hydrophobic interactions', 'Find hydrophobic interactions in the docking poses', self.OnHPhob),
             (wx.ID_ANY, 'Ionic bonds', 'Find ionic bonds in the docking poses', self.OnIonic),
             (wx.ID_ANY, 'Pi-cation', 'Find Pi-cation interactions in the docking poses', self.OnPiCat),
             (wx.ID_ANY, 'Pi-Pi', 'Find Pi-Pi interactions in the docking poses', self.OnPiPi),
             (wx.ID_ANY, 'Sort rescored poses', 'Sort rescored poses', self.OnSortPoses),
             (wx.ID_ANY, 'Virtual Screening - Ligand preparation', 'Prepare ligands for VS', self.OnLigPrep),
             (wx.ID_ANY, 'Virtual Screening - Rename ligands', 'Rename ligands for VS', self.OnLigVS),
             (wx.ID_ANY, 'YOb from mol2', 'Convert mol2 to YOb', self.OnMolYob),
             (wx.ID_ANY, 'YOb to OPM PDB', 'Convert YOb to OPM PDB', self.OnYobOPM),
	      ]:
            if id == None:
                yasaraMenu.AppendSeparator()
            else:
                item = yasaraMenu.Append(id, label, helpText)
                self.Bind(wx.EVT_MENU, handler, item)

        resMenu = wx.Menu()
        for id, label, helpText, handler in \
            [(wx.ID_ANY, 'Plot results with Grace', 'Plot results with grace', self.OnTaskPlot),
             (wx.ID_ANY, 'View Grace graph', 'Open graph in grace', self.OnShowGraph),
             (wx.ID_ANY, 'View results table', 'View results file in viewer', self.OnShowRes),
             (wx.ID_ANY, 'XYCSym', 'Heat plot as XYCSym', self.OnTaskColorbat),
	      ]:
            if id == None:
                resMenu.AppendSeparator()
            else:
                item = resMenu.Append(id, label, helpText)
                self.Bind(wx.EVT_MENU, handler, item)

        viewMenu = wx.Menu()
        for id, label, helpText, handler in \
            [(wx.ID_ANY, 'CMS in PyMOL', 'Show CMS in PyMol', self.OnTaskPyMolCMS),
             (wx.ID_ANY, 'MOL2 in PyMOL', 'Show MOL2 in PyMol', self.OnTaskPyMolMol),
             (wx.ID_ANY, 'PDB in PyMOL', 'Show PDB in PyMol', self.OnTaskPyMolPDB),
             (wx.ID_ANY, 'PDB Sequence', 'Show PDB molecules and chain sequences', self.OnTaskSeq),
             (wx.ID_ANY, 'YOb HBonds', 'List ligand-receptor HBonds in YOb file', self.OnYobHB),
             (wx.ID_ANY, 'YOb Molecules', 'List molecules in YOb file', self.OnYobMol),
	      ]:
            if id == None:
                viewMenu.AppendSeparator()
            else:
                item = viewMenu.Append(id, label, helpText)
                self.Bind(wx.EVT_MENU, handler, item)
                
        appMenu = wx.Menu()
        for id, label, helpText, handler in \
            [(wx.ID_ANY, 'IPython', 'Open IPython console', self.OnTaskIP),
            (wx.ID_ANY, 'Plotting (Grace)', 'Open Grace', self.OnTaskGG),
            (wx.ID_ANY, 'Maestro', 'Open Maestro', self.OnTaskMaestro),
            (wx.ID_ANY, 'VMD', 'Open VMD', self.OnTaskVMD),
            (wx.ID_ANY, 'YASARA', 'Open YASARA', self.OnTaskYASARA),
            (wx.ID_ANY, 'Terminal (Shell)', 'Open shell console', self.OnTaskSh),
	      ]:
            if id == None:
                appMenu.AppendSeparator()
            else:
                item = appMenu.Append(id, label, helpText)
                self.Bind(wx.EVT_MENU, handler, item)

        helpMenu = wx.Menu()
        for id, label, helpText, handler in \
            [(wx.ID_HELP, '&Help', 'Help', self.OnHelp),
             (wx.ID_HELP_INDEX, 'Examples', 'Examples', self.OnExample),
             (wx.ID_ABOUT, '&About', 'Information about this program', self.OnAbout),]:
            if id == None:
                helpMenu.AppendSeparator()
            else:
                item = helpMenu.Append(id, label, helpText)
                self.Bind(wx.EVT_MENU, handler, item)
                
        menuBar = wx.MenuBar()
        menuBar.Append(fileMenu, 'File')
        menuBar.Append(autoDockMenu, 'AutoDock')
        menuBar.Append(DesmondMenu, 'Desmond')
        menuBar.Append(GromacsMenu, 'Gromacs')
        menuBar.Append(pyMenu, 'Python')
        menuBar.Append(vmdMenu, 'VMD')
        menuBar.Append(yasaraMenu, 'YASARA')
        menuBar.Append(resMenu, 'Results')
        menuBar.Append(viewMenu, 'View')
        menuBar.Append(appMenu, 'Applications')
        menuBar.Append(helpMenu, 'Help')
        self.SetMenuBar(menuBar)

    def SetTitle(self):
        # MainWindow.SetTitle overrides wx.Frame.SetTitle, so we have to
        # call it using super:
        super(MainWindow, self).SetTitle('Morelia')

    # Helper method(s):

    def __log(self, message):
        ''' Private method to append a string to the logger text
            control. '''
        self.logger.AppendText('%s\n'%message)

    # Event handlers:
    def OnAbout(self, event):
        dialog = wx.MessageDialog(self, 'Little molecular modelling data analysis app written in python by Jan Jakubik. Non-commercial use allowed under Creative commons licence (BY-NC).', 'About This App', wx.OK)
        dialog.ShowModal()
        dialog.Destroy()

    def OnCWD(self, event):
        dlg = wx.DirDialog(self, "Choose a directory:",
                           style=wx.DD_DEFAULT_STYLE
                           | wx.DD_DIR_MUST_EXIST
                           | wx.DD_CHANGE_DIR
                           )
        if dlg.ShowModal() == wx.ID_OK:
            new_dir = dlg.GetPath()
            self.__log('Current dir is '+str(new_dir))
            os.chdir(new_dir)
            if os.path.isfile("nc.tmp"):
                os.remove("nc.tmp")
            dlg.Destroy()

    def OnPrefs(self, event):
        cmd = ('python '+str(install_dir)+'prefs.py')
        subprocess.Popen(cmd, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        
    def OnHelp(self, event):
        cmd = ('python '+str(install_dir)+'help.py')
        subprocess.Popen(cmd, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)

    def OnExample(self, event):
        dlg = wx.DirDialog(self, "Choose a directory with examples:",
                           defaultPath=example_dir,
                           style=wx.DD_DEFAULT_STYLE
                           | wx.DD_DIR_MUST_EXIST
                           | wx.DD_CHANGE_DIR
                           )
        if dlg.ShowModal() == wx.ID_OK:
            new_dir = dlg.GetPath()
            self.__log('Current dir is '+str(new_dir))
            os.chdir(new_dir)
            dlg.Destroy()

    def OnExit(self, event):
        self.Close()  # Close the main window.

    def OnSave(self, event):
        textfile = open(os.path.join(self.dirname, self.filename), 'w')
        textfile.write(self.logger.GetValue())
        textfile.close()

    def OnEdit(self, event):
        cmd = ('python '+str(install_dir)+'editor.py')
        subprocess.Popen(cmd, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)

    def OnSmi2Pdbqt(self, event):
        self.__log('Converting SMILES to PDBQT ...')
        subprocess.Popen([sys.executable, os.path.join(str(install_dir), 'smi2pdbqt.py')])

    def OnPdb2Pdbqt(self, event):
        self.__log('Preparing receptor PDB to PDBQT ...')
        subprocess.Popen([sys.executable, os.path.join(str(install_dir), 'pdb2pdbqt.py')])

    def OnMapReceptor(self, event):
        self.__log('Mapping receptor with AutoGrid ...')
        subprocess.Popen([sys.executable, os.path.join(str(install_dir), 'map_receptor.py')])

    def OnBatchDock(self, event):
        self.__log('Batch docking with AutoDock-GPU ...')
        subprocess.Popen([sys.executable, os.path.join(str(install_dir), 'batch_dock.py')])

    def OnShowRes(self, event):
        cmd = ('python '+str(install_dir)+'viewer.py')
        subprocess.Popen(cmd, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)

    def OnShowGraph(self, event):
        dialog = wx.FileDialog(self, "Choose Grace file to open", os.getcwd(), "","*.agr", wx.FD_OPEN)
        if dialog.ShowModal() == wx.ID_OK:
            selected=dialog.GetPath()
            cmd = (str(pref.cmdGG)+" "+str(selected))
            subprocess.Popen(cmd, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        dialog.Destroy()
        
    def OnTaskPlot(self, event):
        dialog = wx.FileDialog(self, "Choose a draw file to plot with Grace", os.getcwd(), "","*.draw", wx.FD_OPEN)
        if dialog.ShowModal() == wx.ID_OK:
            selected=dialog.GetPath()
            cmd = (str(pref.cmdGB)+" "+str(selected))
            subprocess.Popen(cmd, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        dialog.Destroy()

    def OnTaskColorbat(self, event):
        self.__log('Heat plot view ...')
        exec_full(str(install_dir)+'colorbat.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnTrjAlign(self, event):
        self.__log('Aligning trajectory ... Please wait ...')
        exec_full(str(install_dir)+'Trj_Align.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnLigRMSD(self, event):
        self.__log('Ligand RMSD ...')
        exec_full(str(install_dir)+'Lig_RMSD.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnPLCon(self, event):
        self.__log('Protein-Ligand Contacts ...')
        exec_full(str(install_dir)+'PL_Contacts.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnMMin(self, event):
        subprocess.Popen('$SCHRODINGER/run desmond_minimization_gui.py', shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)

    def OnMD(self, event):
        subprocess.Popen('$SCHRODINGER/run desmond_molecular_dynamics_gui.py', shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)

    def OnDummy(self, event):
        self.__log('Dummy GPU-host to local host ...')
        exec_full(str(install_dir)+'dummy.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnPBBRMSD(self, event):
        self.__log('Protein backbone RMSD ...')
        exec_full(str(install_dir)+'Prot_BB_RMSD.py')

    def OnPRMSF(self, event):
        self.__log('Protein RMSF ...')
        exec_full(str(install_dir)+'Prot_RMSF.py')

    def OnSID(self, event):
        subprocess.Popen('$SCHRODINGER/run event_analysis.py', shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)

    def OnQA(self, event):
        subprocess.Popen('$SCHRODINGER/run simulation_analysis_gui.py', shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)

    def OnTrjRMSD(self, event):
        subprocess.Popen('$SCHRODINGER/run trajectory_cluster_gui.py', shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)

    def OnLigDistance(self, event):
        self.__log('Measuring ...')
        exec_full(str(install_dir)+'measure_ligand_centroid_wx_app.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnTRR2XTC(self, event):
        self.__log('Compressing .trr trajectory to .xtc ...')
        subprocess.Popen([sys.executable, os.path.join(str(install_dir), 'trr2xtc.py')])

    def OnDesmondMMPBSA(self, event):
        self.__log('Starting Desmond MM-PBSA ...')
        subprocess.Popen([
            sys.executable,
            os.path.join(
                str(install_dir),
                'desmond_gmxMMPBSA',
                'desmond_gmxmmpbsa_wx.py',
            ),
        ])

    def OnCMS2GRO(self, event):
        self.__log('Converting .cms structure to .gro ...')
        subprocess.Popen([
            sys.executable,
            os.path.join(str(install_dir), 'cms2gro.py'),
        ])

    def OnDTR2TRR(self, event):
        self.__log('Converting .dtr trajectory to .trr ...')
        subprocess.Popen([
            sys.executable,
            os.path.join(str(install_dir), 'dtr2trr.py'),
        ])

    def OnArchive(self, event):
        self.__log('Achiving ...')
        exec_full(str(install_dir)+'archive.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnDockLS(self, event):
        self.__log('Docking ...')
        exec_full(str(install_dir)+'dock_local.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnDockRe(self, event):
        self.__log('Rescoring ...')
        exec_full(str(install_dir)+'dock_rescore.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnExtractTop(self, event):
        self.__log('Extracting top poses ...')
        exec_full(str(install_dir)+'extract_top_poses.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnInteractions(self, event):
        self.__log('Analyzing interactions ...')
        exec_full(str(install_dir)+'interactions.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnMakeSiTopo(self, event):
        self.__log('Making Si ligand topology ...')
        dialog = wx.FileDialog(self, "Choose Si ligand file", os.getcwd(), "", "SDF/MOL2 files (*.sdf;*.mol2)|*.sdf;*.mol2", wx.FD_OPEN)
        if dialog.ShowModal() == wx.ID_OK:
            selected=dialog.GetPath()
            subprocess.Popen([sys.executable, os.path.join(str(install_dir), 'make_usertopo.py'), selected])
        dialog.Destroy()

    def OnMinSiLigands(self, event):
        self.__log('Minimizing Si ligands ...')
        subprocess.Popen([sys.executable, os.path.join(str(install_dir), 'minimize_si_ligands.py')])

    def OnHBonds(self, event):
        self.__log('Hydrogen bonds ...')
        exec_full(str(install_dir)+'H_Bonds.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnIonic(self, event):
        self.__log('Ionic bonds ...')
        exec_full(str(install_dir)+'ionic.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnHPhob(self, event):
        self.__log('Hydrophobic interactions ...')
        exec_full(str(install_dir)+'HPhob.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnPiCat(self, event):
        self.__log('Pi-cation interactions ...')
        exec_full(str(install_dir)+'PiCat.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnPiPi(self, event):
        self.__log('Pi-Pi interactions ...')
        exec_full(str(install_dir)+'PiPi.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnSortPoses(self, event):
        self.__log('Sorting of poses ...')
        exec_full(str(install_dir)+'sort_poses.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnLigPrep(self, event):
        self.__log('Preparing ligands ...')
        exec_full(str(install_dir)+'mol2_em_yob.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnLigVS(self, event):
        self.__log('Renaming ligands ...')
        exec_full(str(install_dir)+'vs_rename.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnMolYob(self, event):
        self.__log('Converting ligand ...')
        exec_full(str(install_dir)+'mol2_em_yob_single.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnYobOPM(self, event):
        self.__log('Converting ...')
        exec_full(str(install_dir)+'convert_yob_to_OPM_pdb.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnPDBOPM(self, event):
        self.__log('Converting ...')
        exec_full(str(install_dir)+'align_pdb_to_membrane.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnTaskPyMolCMS(self, event):
        self.__log('Showing CMS 3D structure ...')
        exec_full(str(install_dir)+'cms_pymol.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnTaskPyMolMol(self, event):
        self.__log('Showing ligand 3D structure ...')
        exec_full(str(install_dir)+'mol2_pymol.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnTaskPyMolPDB(self, event):
        self.__log('Showing PDB 3D structure ...')
        exec_full(str(install_dir)+'pdb_pymol.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)

    def OnTaskSeq(self, event):
        self.__log('Reading PDB ...')
        exec_full(str(install_dir)+'pdb_seq.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)
        log = open('pdb_seq.res','r')
        tolog = log.read()
        self.__log(tolog)

    def OnYobMol(self, event):
        self.__log('Reading YOb ...')
        exec_full(str(install_dir)+'list_mol_in_yob.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)
        log = open('mol_list.res','r')
        tolog = log.read()
        self.__log(tolog)

    def OnYobHB(self, event):
        self.__log('Reading YOb ...')
        exec_full(str(install_dir)+'list_HBonds_in_yob.py')
        log = open('temp.log','r')
        tolog = log.read()
        self.__log(tolog)
        log = open('HBonds.tab','r')
        tolog = log.read()
        self.__log(tolog)

    def OnTaskIP(self, event):
        subprocess.Popen(pref.cmdIP, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)

    def OnTaskGG(self, event):
        subprocess.Popen(pref.cmdGG, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)

    def OnTaskYASARA(self, event):
        subprocess.Popen(pref.cmdYASARA, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)

    def OnTaskSh(self, event):
        subprocess.Popen(pref.cmdSh, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)

    def OnTaskMaestro(self, event):
        subprocess.Popen(pref.cmdMaestro, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)

    def OnTaskVMD(self, event):
        subprocess.Popen(pref.cmdVMD, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)

    def OnTask(self, event):
        dialog = wx.MessageDialog(self, 'Not implemented yet.', 'Ehm', wx.OK)
        dialog.ShowModal()
        dialog.Destroy()
        

class DataSets(wx.Dialog):
    def __init__(
            self, parent, ID, title, size=wx.DefaultSize, pos=wx.DefaultPosition, 
            style=wx.DEFAULT_DIALOG_STYLE,
            useMetal=False,
            ):
        wx.Dialog.__init__(self)
        self.SetExtraStyle(wx.DIALOG_EX_CONTEXTHELP)
        self.Create(parent, ID, title, pos, size, style)

        sizer = wx.BoxSizer(wx.VERTICAL)
        label = wx.StaticText(self, -1, "How many data sets do you want to analyse?")
        sizer.Add(label, 0, wx.ALIGN_CENTRE|wx.ALL, 5)
        
        sampleList = ['1', '2', '3', '4', '5']
        rb = wx.RadioBox(
                self, -1, "", wx.DefaultPosition, wx.DefaultSize,
                sampleList, 5, wx.RA_SPECIFY_COLS | wx.NO_BORDER
                )
        self.Bind(wx.EVT_RADIOBOX, self.EvtRadioBox, rb)
        rb.SetToolTip(wx.ToolTip("Select the number"))
        sizer.Add(rb, 0, wx.LEFT|wx.RIGHT|wx.BOTTOM, 20)

        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
        
        btnsizer = wx.StdDialogButtonSizer()
        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()
        btnsizer.AddButton(btn)
        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_RIGHT|wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)
    def EvtRadioBox(self, event):
        choice = (event.GetSelection())
        curves = float(choice) + 1
        nc = open('nc.tmp','w')
        nc.write(str(curves))
        nc.close()
       
app = wx.App(False)
frame = MainWindow()
frame.Show()
app.MainLoop()
