#!/opt/anaconda-2.0/bin/python
"""
The jobsubmit.py modules provides methods for submitting and running
jobs that use qsub and PBS directives.

jobWait (boolean): If true. when trying to send more than total number
    of jobs
"""
#python 2.x and 3.x cross compatibility
from __future__ import print_function, division
import sys
if sys.version_info[0] >= 3:
    xrange = range

import os
import time
from subprocess import Popen, PIPE, call


jobWait =False #set to false if running on webserver

def IsUserOverJobLimit(username,maxJobs):
    """
    Checks if user over maximum number of jobs
    Arguments:
        username (str): name of current user
        maxJobs (int): Maximum number of jobs allowed to be submitted
    Returns:
        True if over maximum limit. False otherwise.
    """
    maxJobs = int(maxJobs)
    p = Popen("qstat |grep "+username + " |wc -l",shell = True,
                stdout=PIPE,stderr=PIPE)
    stdout, stderr = p.communicate()
    if int(stdout) > maxJobs:
        return True
    else:
        return False


class JOBS:
    """
    This class is for writing and sending jobs
    Attributes:
        username (str): username of the user submitting jobs
        walltime (str): hrs:min:sec Maximum length for job to run
        priority (str): Submission priority. choices are: default and urgent
        maxJobs (int): Maximum number of jobs allowed to be submitted at once.
    Methods:
        SubmitJob: Submits a Job.
        CreateJob: Creates a shell script for a job
        JobIsDoneOrRunning: Checks if job is completed or running
        Execute: Runs an executable command locally. Prints stdout in realtime.
            Returns stdout and stderror
        mkdir: Attempts to make a given directory
        KillJob: Kill a particular job
        KillAll: Kill all jobs of a user, user needs to gointo function
            and delete protection that ignores users except USER
        CleanTmp: Clean all temporary directories of a user on all nodes
    """

    def __init__(self,username,walltime="24:59:00", priority='default',
                    maxJobs = 300 ):
        self.username = username
        self.walltime = walltime
        self.priority = priority
        self.maxJobs = maxJobs

    def SubmitJob(self,jobInput, jobOutput, jobName, outputDir, PBSerror = '',
        PBSoutput = '',ForceSubmit = False, WaitIfOverJobLimit=jobWait):
        """
        Creates and submits a job if the job is not: in que, running or
        completed.

        Arguments:
            jobInput ([str]): List of commands use for running the script.
                ie. jobInput = ['/path/exampleJob', 'arg1','arg2'].
                Each item in the list should contain the full path. 
            jobOutput ([str]): Full path name for all files and directories
                that should be returned after job is completed. Used
                to for checking if job has finished.
            jobName (str): Name of job to be submitted.
            outputDir (str): Full path of output directory
            PBSerror (str,optional): Full path name of pbs error file. Default
                uses outputDir/jobName.err
            PBSoutput (str,optional): Full path name of pbs output file. Default
                uses outputDir/jobName.out
            ForceSubmit (bool, optional): Default is False.  Used to allow
                user to force a job submission if the output files 
                are present.  The output may be present but not complete
                if user detects this can force resubmission.
            WaitIfOverJobLimit: Puts system to sleep in over job limit. If false
                skips submitting job if over limit.
        """
        if self.JobIsRunning(jobName):
            return


        if not self.JobIsDone(jobName, jobOutput) or ForceSubmit:
            self.CreateJob(jobInput, outputDir, jobName, PBSerror = '',
                PBSoutput = '')
            #check number of jobs running, if over max wait for some jobs to finish
            while True:
                p = Popen("qstat |grep "+self.username + " |wc -l",shell = True,
                            stdout=PIPE,stderr=PIPE)
                stdout, stderr = p.communicate()
                if int(stdout) > self.maxJobs:
                    if WaitIfOverJobLimit:
                        time.sleep(60)
                    else:
                        print(jobName,'submission skipped due to being over job limit')
                        return
                else:
                    break    
            #run job
            p = Popen('qsub '+outputDir + '/' +jobName, shell=True,stdout=PIPE,
                            stderr=PIPE)
            stdout, stderr = p.communicate()
            out = stdout.split('\n')
            for line in out:
                if line != '':
                    print(line.strip())
            print(jobName, "submitted")
            time.sleep(0.1)

    def CreateJob(self,jobInput, outputDir, jobName, PBSerror = '',
                    PBSoutput = ''):
        """
        Creates a shell PBS script for job submission.

        Arguments:
            jobInput ([str]): A list containing the inputs for the job
                script. ie jobInput = ['./exampleJob.py','arg1','arg2',etc.]
            outputDir (str): The full path for the output directory.  Unless
                otherwise specified the job script, output and error files
                are written there.
            jobName (str): The name of the job without file extensions. ie
                'exampleJob' is okay but 'exampleJob.txt' not okay.
            PBSerror (str,optional) Full path name of PBS error file.  Default
                is the jobName.err stored in the outputDir
            PBSoutput (str,optional): Full path name of PBS output file. Default
                is the jobName.out stored in the outputDir 
        Returns:
            directory (directory): Attempts to create the output directory.
            PBS script (file): PBS job submission script. Which is written to
                outputDir/jobName   
        """
        #job inputs and outputDir require full paths
        if PBSerror == '':
            PBSerror = self.GetPbsPath(outputDir,jobName)[0] #outputDir + '/' + jobName + ".err"
        if PBSoutput =='':
            PBSoutput =self.GetPbsPath(outputDir,jobName)[1] #outputDir + '/' + jobName + ".out"
        self.mkdir(outputDir)
        fjob = open(self.GetPbsPath(outputDir,jobName)[2], 'w')
        fjob.write("#PBS -q "+self.priority +'\n')
        fjob.write("#PBS -d /tmp/"+'\n')
        fjob.write("#PBS -l nodes=1:ppn=1,walltime="+self.walltime+'\n')
        fjob.write("#PBS -o "+PBSoutput +'\n')
        fjob.write("#PBS -e "+PBSerror + '\n')
        fjob.write("#PBS -N "+jobName +'\n')    
        for input in jobInput:
            fjob.write(str(input))
            fjob.write(' ')
        fjob.write('\n')
        fjob.close()

    def GetPbsPath(self,outputDir,jobName):
        """
        Returns PBS output and error file paths
            PBSerror  (str) Full path to pbs error file
            PBSoutput (str) Full path to pbs output file
            PBSsubmit (str) Full path to pbs submision file
        """
        PBSerror = outputDir + '/' + jobName + ".err"
        PBSoutput= outputDir + '/' + jobName + ".out"
        PBSsubmit = outputDir+'/'+jobName
        return PBSerror, PBSoutput,PBSsubmit

    def JobFailed(self,outputDir,jobName,outputFileList):
        """
        Checks to see if job failed. Returns true if pbs
        output or error file exists but not all outputFiles
        have been generated
        """
        PBSerror,PBSoutput,PBSsubmit = self.GetPbsPath(outputDir,jobName)
        jobFinish = os.path.exists(PBSerror) or os.path.exists(PBSoutput)
        if not jobFinish:
            return False
        for outputFilePath in outputFileList:
            if not os.path.exists(outputFilePath):
                return True
        return False

    def JobIsDone(self,jobName, jobOutput = None):
        """
        Checks if the job has been completed by checking for output files.

        Arguments:
            jobName (str): Name of job.
            jobOutput ([str]): List of output files and directories from job.  
                The program just checks for existance.  It does not check
                if the files have actually completed. If list is empty
                function returns False.
        Returns:
            True: If all output files are generated
            False: otherwise
        """ 
        if jobOutput is None:
            return False
            
        if not jobOutput:
            return False

        for output in jobOutput:
            if not os.path.exists(output):
                return False
        print(jobName, "completed skipping submission")
        return True


    def JobIsRunning(self,jobName,verbose=True):
        """
        Checks if the job is currently running.
        Arguments:
            jobName (str): Name of job.
            verbose (bool): If true prints statements
        Returns:
            True: If job is running or in que.
            False: otherwise
        """ 
        p = Popen("qstat -f|grep Job_Name", shell = True, stdout = PIPE,
                    stderr = PIPE) 
        stdout, stderr = p.communicate()
        out = stdout.split('\n')
        for name in out:
            if jobName in name:
                if verbose:
                    print(jobName, "is running or in que skip submission")
                return True
        return False


    def Execute(self, command, raiseStdError = True):
        """
        Runs a system command.  Prints stdout in real time. Returns stdout and 
        stderr.
        
        Arguments:
            command (str or [str]): System command can be run.  Can be a string 
                or a list.  ie command = 'ls -la' or command = ['ls','-la']
            raiseStdError (bool): Default True causes an error to be raised if 
                the executable writes output to standard error. If set to false.
                stderr is ignored.
        Returns:
            stdout (str): Standard output written by executable.
            stderr (str): Standard error written by program. 
        """
        if isinstance(command, list):
            cmd = ' '.join(command)
        else:
            cmd = command   
        p = Popen(cmd, shell=True, stdout=PIPE, stderr = PIPE)
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


 

    def mkdir(self, directory):
        try:
            os.mkdir(directory)
        except:
            pass
    def KillJob(self,jobName='',jobID=''):
        """
        Kills the given job.

        Requires:
            jobName or JobID to be given
        Arguments:
            jobName (str,optional): Name of job to be killed.
            jobID   (str or integer,optional): Job ID of job to be killed  
        """
        if jobID != '':
            print(jobID, "killed")
            call(['qdel','-p',str(jobID)])
        elif jobName != '':
            print(jobName, "killed")
            call(['qdel','-p',jobName])
        else:
            pass    
        
    def KillAll(self):
        """
        Kill all jobs for a particular user.
        
        Arguments:
            username (str): The user name of the user. 
        """
        if self.username != os.getenv("USER"):
            strout = self.username +" Not approved for this function you can change KillAll function to allow"
            strout += " user to use KillAll function"
            print(strout)
            return
        p=Popen("qstat -u "+self.username, shell = True,stdout=PIPE,stderr=PIPE)
        stdout, stderr = p.communicate()
        out = stdout.split()
        joblist = []
        for line in out:
            if line == '':
                continue
            if line[0].isdigit():
                jobinfo = line.split()
                jobID = jobinfo[0].split('.')[0]
                jobName = jobinfo[3]
                joblist.append([jobID,jobName])
        for job in joblist:
            self.KillJob(jobID = job[0] )
    def CleanTmp(self):
        """
        Delete all /tmp/username directories and its subdirectories
        on each node.

        Arguments:
            username (str): The user name for the user.
        """
        call(["tentakel", "/bin/rm", "-rf", "/tmp/"+self.username])
