#configuration dictionary for hhsearch
import sys
import os

#file paths and username
hhsearchDirectory = os.path.dirname(os.path.realpath(__file__) ) +'/'
sourceDirectory = os.path.dirname(os.path.realpath(hhsearchDirectory))+'/'
sys.path.append(sourceDirectory)
import userConfig

pythonPath = userConfig.pythonPath
parameters = {
'seqfileExtension':'.seq',                  #extension to query sequence file if file is 117eA.fasta ext is .fasta common ext [.fasta,.seq,'']
'rootDir':sourceDirectory,                  #Path to source root. ie /nfs/amino-home/bgovi/source/
'pdbPath':sourceDirectory+'/PDB/',              #Path to PDB directory
'jobPath':sourceDirectory+'/JOB/',              #Path to job directory
'commandPath':sourceDirectory+'/CommandFunctions/',     #Path to command functions directory
}

