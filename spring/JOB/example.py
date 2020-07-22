#!/opt/anaconda-2.0/bin/python
#simple example script for submitting a job
import os
import jobsubmit

def mkdir(directory):
    try:
        os.mkdirs(directory)
    except:
        pass

#change these variables
username = 'bgovi'
scriptPath = "/nfs/amino-home/bgovi/source/JOB/example/test.py"
jobTestOutput = "/nfs/amino-home/bgovi/source/JOB/example/job_test.txt"
jobName = 'JobTest'
outputDir = "/nfs/amino-home/bgovi/source/JOB/example"
mkdir(outputDir)

job = jobsubmit.JOBS(username)  #initializes module.

"""
The job object requires 4 variables for job submission.
jobInput is the program you want to run and its command line arguments.
    jobInput = ['/pathToprogram/program','argument1','argument2',...]

jobOutput is a list of all the output files.  If all the files are completed
the jobsubmission is skipped.  If list is empty job will be submitted.
    jobOutput = ['/fullPath/ToOutputFile1','/fullPath/ToOutputFile2',..]

jobName is the job name which should be unique.  The module will check to
see if the program is currently submitted or running and if so will skip 
submission.

outputDir is where the pbs script the standard error and standard output files
from the jobsubmission are written to.
"""



jobInput = [scriptPath]
jobOutput = [jobTestOutput]

job.SubmitJob(jobInput, jobOutput, jobName, outputDir)
