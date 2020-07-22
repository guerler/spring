"""
Please change the username variable to your username.

The majority of directories in source use this file
to get the username.  The position of this file within
the source directory can not be changed.

The getpass module should return the current users name
if it fails you can manually change the user name here.


There are two protocols for jobs when trying to submit over the
limit.  You can either wait or kill the script. On the webserver
the script should die. So in the JOB/jobsubmit.py file if running
on webserver switch jobWait variable to False.
"""
import getpass
import os,sys

sourceDirectory = os.path.dirname(os.path.realpath(__file__))+'/'

username = getpass.getuser()    #'ie bgovi,if user is bgovi. '
pythonPath ='/usr/bin/env python'#'/opt/anaconda-2.0/bin/python'
temporaryDirectory = '/tmp/'    #root location of where all temporary files are written
emailAddress=''         #Optional user email
webserver=''

##Path for databases##
#rootdbPath='/nfs/amino-home/bgovi/database/'    #path to where all databases are stored
rootdbPath=sourceDirectory+'database/'    #path to where all databases are stored
#aminoLibrary='/nfs/amino-library/'              #local directory on all nodes
aminoLibrary=rootdbPath+'amino-library/'              #local directory on all nodes
#aminoZhangUniprotDir = aminoLibrary+'local/hhsuite/uniprot20_2015_06/'
aminoZhangUniprotDir = aminoLibrary+'local/hhsuite/uniprot20_2016_02/'
#rootdbPath='/nfs/amino-library/DIMERDBX/'
#Variables for parsing PDB#
substrPos=1
substrLength=2
"""
substrPos and substrLength are variables for parsing the pdb directory format.
For example the pdb chain 12asA in directory pdb will be stored in
./pdb/2a/12asA.pdb.  The 2a comes from the second position and a substring
of length 2.  This file structure is used in most of the databases here.
"""

#Job Parameters Default.  Can be change for each module
#in its respective configuration file
maxJobs = 1200
priority = 'default'
waitIfOverJobLimit=True
"""
If over job limit wait in background or kill script.
Need to set to false when running on webserver
"""

#PDB Database
pdbPath=rootdbPath+'pdb/'       #all pdb data is stored here
pdbAllPath=pdbPath+'PDBall/'        #directory contains all protein complexes. Extract from pdbPath
PDBNonRedundant=pdbPath+'PDBNonRedundant/' #Contains list of non redundant complexes contains in pdbAllPath
chainsPath=pdbPath+'chains/'        #Contains all unique monomers from each pdbid from the pdbPath directory
seqPath=pdbPath+'sequences/'        #Contains fasta format sequences of proteins in chainsPath
pdbDatabase = pdbPath+'pdb/'        #all pdb files downloaded from www.rcsb.org

#Databases for HHsearch
HHsearchdbPath=rootdbPath+'HHsearch/'   #Path to where the hhsearch databases are stored
hhmDB=aminoLibrary+'/DIMERDB/HHsearch/hhm.db'       #Path to hhm profile database
uniprotDir=aminoZhangUniprotDir  #Path to uniprot database
#uniprotDBName='uniprot20_2015_06'   #Name of uniprot database
uniprotDBName='uniprot20_2016_02'   #Name of uniprot database

#Data profiles loctation for COTH
cothdbPath=rootdbPath+'COTH/'       #Path to the COTH profiles

#SPRING index files locations
springDB=rootdbPath+'SPRING/'       #Path to the SPRING index files
springIndexAll=springDB+'indexAll.txt'  #List for spring contains all pairwise protein complexes
springIndex=springDB+'index.txt'    #Non redundant list from indexAll.txt

#Blast database location
blastdbPath=rootdbPath+'BLAST/'     #Path to blast databse directory
blastdb=blastdbPath+'blastPDB.seq'  #Full path to blast database
if __name__ == "__main__":
    print username
