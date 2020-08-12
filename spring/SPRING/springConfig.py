import sys
import os

springDirectory = os.path.dirname(os.path.realpath(__file__) ) +'/'
sourceDirectory = os.path.dirname(os.path.realpath(springDirectory))+'/'
sys.path.append(sourceDirectory)
spring_db = '/home/bgovi/Workspace/SPRING/spring/spring/spring.test-data/springdb/'
tmp_dir = '/home/bgovi/Workspace/SPRING/spring/spring/tmp/'


#configuration dictionary for spring.py
parameters = {
'user': 'bgovi',
'tmpDir': tmp_dir,
'rootDir':sourceDirectory,      #Path to source root. e.g. $HOME/source/
'springPath':springDirectory,       #Full path to SPRING directory
'pdbPath':sourceDirectory+'/PDB/',      #Path to PDB subdirectory
'commandPath':sourceDirectory+'/CommandFunctions/',     #Path to CommandFunction subdirectory
'hhsearchPath':sourceDirectory+'/HHsearch/',            #Path to HHsearch Directory
'maxNumberTemplates':100,       #Maximum Number of Tempaltes returned by SPRING
'maxNumberTemplateInit':100,       #Maximnum Number of Templates in init.dat. This file used by TACOS simulation
'WriteFullSummary':True,        #write all templates to summary even if over max templates
'seqidCutoff':1.1,#1.1,         #sequence id cutoff.  benchmark usually set at 0.3.Real mode above 1.0
'seqFileExtension':'.seq',      #.fasta or .seq. Depending on sequnce extension used in target directories
"W0": 12.0,                 #Weight for SPRING Score
"W1": 1.4,                  #Weight for SPRING Score
"Normalize": 6.5,           #Weight for SPRING Score
"runHHsearch": False,           #run HHsearch if hhr files missing.
'indexFile': spring_db+'index.txt', #Full path to index.txt file
'dfireFile':springDirectory+'/dfire.txt',   #Full path to dfire.txt file
'homodimerThreshold': 0.70,  #Homodimer threshold. If targetA aligned to targetB have seqid below threshold run COTH in both orders
'maxModels':50,
}

substrPos = 0
substrLength = 2

PDBparameters = {
'pdbDir': spring_db +'pdball/',      #Full Path to root of PDB complex Structures
'monomerpdbDir':spring_db +'chains/',        #Full path to all monomer PDBs that are in the HHsearch database
'monomerPDBExtension': '.pdb',  #Normally .pdb, but the current PDB database for hhsearch does not have the .pdb extension so left empty
'substrPosition': substrPos,
'substrLength': substrLength,
'substrEnd':substrPos + substrLength,
'seqFileExtension':'.seq',  #Target sequence extension for sequences written in fasta format
'pdbFileExtension':'.pdb',  #PDB extension for pdbs in the complex library
}

