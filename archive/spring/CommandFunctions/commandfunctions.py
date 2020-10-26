#!/usr/bin/env python
#python 2.x and 3.x cross compatibility
from __future__ import print_function, division
import sys
if sys.version_info[0] >= 3:
    xrange = range

import os
import shutil
from subprocess import Popen, PIPE

def Execute(command, raiseStdError = True,environment = None,verbose=True):
    """
    Runs a system command.  Prints stdout in real time. Returns stdout and 
    stderr.
    
    Arguments:
        command (str or [str]): System command can be run.  Can be a string or 
            a list.  ie command = 'ls -la' or command = ['ls','-la']
        raiseStdError (bool): Default True causes an error to be raised if the 
            executable writes output to standard error.  If set to false.
            stderr is ignored.
        enviornment (dictionary): Default enviornment is os.environ
        verbose (bool): Print output to screen in real time.
    Returns:
        stdout (str): Standard output written by executable.
        stderr (str): Standard error written by program. 
    """
    if isinstance(command, list):
        cmd = ' '.join(command)
    else:
        cmd = command   

    if environment is None:
        myenv = os.environ
    else:
        tmpEnv = os.environ.copy()
        tmpEnv.update(environment)
        myenv = tmpEnv

    if verbose:
        p = Popen(cmd, shell=True, stdout=PIPE, stderr = PIPE,env=myenv)
        stdoutTmp = []

        for line in iter(p.stdout.readline, ''):
            print(line,end='')
            stdoutTmp.append(line)
        stderTmp = []
        for line in p.stderr:
            stderTmp.append(line)
        stdout = ''.join(stdoutTmp)
        stderr = ''.join(stderTmp)
        if raiseStdError and stderr != '':
            print(command, "Prints error messages")
            print(stderr)
            raise RuntimeError
        elif stderr != '':
            print(stderr)
        return stdout, stderr
    else:
        p = Popen(cmd, shell=True, stdout=PIPE, stderr = PIPE,env=myenv)
        stdout, stderr = p.communicate()
        if raiseStdError and stderr != '':
            print(command, "Prints error messages")
            print(stderr)
            raise RuntimeError
        elif stderr != '':
            print(stderr)
        return stdout, stderr
def mkdir(directory):
    """
    Makes directory.  Ignores error such as path already existing
    """
    try:
        os.mkdir(directory)
    except:
        pass

def mkdirs(directory):
    """
    Makes all intermediate directories needed to contain the leaf directory
    """
    try:
        os.makedirs(directory)
    except:
        pass
def RemoveDirectory(directory):
    """
    Deletes directory and all subdirectories.  Does not allow
    deletion of root directory.

    Arguments:
        directory (str): Full path of directory to be deleted
    """
    directory = directory.strip()
    #check if directory is root directory ie '/'
    if '' == directory:
        return
    if '' == ''.join(directory.split('/')):
        print("Tried to remove root directory, skipping removal of directory")
        return
    else:
        shutil.rmtree(directory)            
