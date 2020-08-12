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

if __name__ == "__main__":
    print(username)
