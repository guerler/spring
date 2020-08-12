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
'user':userConfig.username,                     #username ie. bgovi
'tmpDir':userConfig.temporaryDirectory,
'seqfileExtension':'.seq',                  #extension to query sequence file if file is 117eA.fasta ext is .fasta common ext [.fasta,.seq,'']
'rootDir':sourceDirectory,                  #Path to source root. ie /nfs/amino-home/bgovi/source/
'hhsearchPATH':hhsearchDirectory,               #Path where hhsearch.py is located
'pdbPath':sourceDirectory+'/PDB/',              #Path to PDB directory
'jobPath':sourceDirectory+'/JOB/',              #Path to job directory
'commandPath':sourceDirectory+'/CommandFunctions/',     #Path to command functions directory
'hhsuiteDir':hhsearchDirectory+"/hhsuite/",             #Path location of hhsuite directory.  Contains HHsearch program
'uniprotDir':userConfig.uniprotDir,                 #Database location for hhmake
'uniprotDBName':userConfig.uniprotDBName,           #uniprot DB name
'hhsearchDB':userConfig.hhmDB,                  #Database location for hhsearch
'hhblits':hhsearchDirectory+"/hhsuite/bin/hhblits",         #hhblits Path
'hhmake':hhsearchDirectory+"/hhsuite/bin/hhmake",       #hhmake path
'hhsearch':hhsearchDirectory+"/hhsuite/bin/hhsearch",       #hhsearch path
'addss.pl':hhsearchDirectory+"/hhsuite/scripts/addss.pl",   #path to scripts where addss.pl is located
'runMake':False,                        #run hhmake [True,False]
'runSearch':False,                      #run hhsearch [True,False]
'runComplex':False,                     #run complex ie 117eA or 117eA-117eB [True,False]
'runMonomer':False,                     #for running monomers ie 117eA.  Where each monomer has individual target directory
'runGenome': False,                     #for running monomers with the directories in genome wide format.
#GenomeWide Prediction parameters for directories
'substrPosition':userConfig.substrPos,
'substrLength':userConfig.substrLength
}

hhblitsParameters = {
'-d':parameters['uniprotDir'] + '/' + parameters['uniprotDBName'], #database path and database name
'-n':'2',                           #number of interations [1,8]
'-e':'0.001'                            #E-value cutoff for inclusion in result alignment [0,1]
}

hhmakeParameters = {
'-id': '90', #maximum pairwise sequence identity (%) [0,100]
'-diff': '100', #filter most diverse set of sequences, keeping at least this many sequence in each block of >50 columns
'-cov': '0', #minimum coverage with query (%) [0,100]
'-qid': '0' #minimum sequence identity with query (%) [0,100]
}

hhsearchParameters = {
'-d': parameters['hhsearchDB'], #hhsearch database
#Arguments for filtering
'-id': '90',    #maximum pairwise sequence identity (%) [0,100]
'-diff': '100', #filter most diverse set of sequences, keeping at leaasst this many sequence in each block of >50 columns
'-cov': '0',    #minimum coverage with query (%) [0,100]
'-qid': '0',    #minimum sequence identity with query (%) [0,100]
#Arguments for output
'-e':'0.001',   #E-Value cutoff for inclusion in multiple alignment
'-p':'20',  #minimum probability in summary and alignment
'-E':'0.01',    #maximum E-Value in summary and alignment
'-Z':'30000', #maximum number of lines in summary hit
'-z':'20000',  #minimum number of lines in summary hit
'-B':'30000',   #maximum number of alignments in alignment list
'-b':'20000'   #minimum number of alignments in alignment list
}

jobParameters = {
'batchSize':1,    #number of targets run per job
'walltime':"50:59:00",  #Job time limit "hours:minutes:seconds"
'priority':userConfig.priority,#Job submission priority.  ie default,urgent,casp
'maxJobs': userConfig.maxJobs,     #maximum number of jobs to submit
'forceSubmit':False #force resubmission when all results files have already been generated.
}
