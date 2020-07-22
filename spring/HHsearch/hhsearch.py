#!/usr/bin/env python
"""
hhsearch.py is a program for running and submitting hhsearch programs.
hhsearchConfig.py is a configuration file with paths and parameters that may
need to be altered.
"""
#python 2.x and 3.x cross compatibility
from __future__ import print_function, division
import sys
if sys.version_info[0] >= 3:
    xrange = range

import os
import sys
import os.path
import shutil
import argparse

import hhsearchConfig
rootDir = hhsearchConfig.parameters['rootDir']
sys.path.append(hhsearchConfig.parameters['jobPath'])
from jobsubmit import JOBS
sys.path.append(hhsearchConfig.parameters['commandPath'])
import commandfunctions

#Need to create environment for running hhsearch system commands.
os.environ["HHLIB"] = hhsearchConfig.parameters['hhsuiteDir']

class HHSEARCH:
    """
    This class is a module for wrapping the template threading program hhsearch.
    The module provides wrappers for hhmake and hhsearch.  

    Attributes (default parameters in hhsearchConfig.py):
        query: Full path of newline delimited query targets or name of query 
               target
        inputDirectory: Directory containing fasta,hhm and hhr directory for 
                single chain run, or directory where all complex
                target directories are. ie if directories 117eA-117eB,
                12asA-12asB are in /home/bgovi/test/ then 
                inputDirectory = /home/bgovi/test/
        user: Name of user
        hhsearchPATH: Full path to directory where hhsearch.py is located.
        hhsuiteDir: Full path to hhsuite directory
        uniprotDir: Full path to location of uniprot database
        uniprotDBName: Name of uniprot database
        hhserachDB: Full path of hhsearch database.
        hhblits: Full path for hhblits executable
        hhmake: Full path for hhmake executable
        addss: Full path for addss.pl perl script
        runMake: Boolean for running hhmake
        runSearch: Boolean for running hhsearch
        seqfileExtension: Sequence files have several default extensions suchas 
                '.fasta','.seq',''.  Used to choose which extension used for 
                query sequence
        tmpDir: temporary directory for running hhmake/hhsearch.

    Methods:
        HHsearchMain:   Main program 
        GetAndMakePaths: Given a query and inputDirectory returns locations of 
                fasta,hhm,and hhr directories. (the directory is different for 
                complexes ie 117eA-117eB than for 117eA and 117eB seperately.
                read README for details.)
        RunHHmake_HHsearch: Runs HHmake and HHsearch
        Make:   Wrapper function for running hhmake 
        Search: Wrapper function for running hhsearch
        Submit: For submitting hhsearch and/or hhmake jobs

    """
    def __init__(self,query,inputDirectory,
            user =hhsearchConfig.parameters['user'], 
            hhsearchPATH = hhsearchConfig.parameters['hhsearchPATH'],
            hhsuiteDir = hhsearchConfig.parameters['hhsuiteDir'], 
            uniprotDir=hhsearchConfig.parameters['uniprotDir'],
            uniprotDBName = hhsearchConfig.parameters['uniprotDBName'], 
            hhsearchDB = hhsearchConfig.parameters['hhsearchDB'],
            runMake = hhsearchConfig.parameters['runMake'], 
            runSearch = hhsearchConfig.parameters['runSearch'],
            seqfileExtension = hhsearchConfig.parameters['seqfileExtension'],
            runComplex = hhsearchConfig.parameters['runComplex'],
            runMonomer = hhsearchConfig.parameters['runMonomer'], 
            runGenome = hhsearchConfig.parameters['runGenome'],
            substrPosition = hhsearchConfig.parameters['substrPosition'], 
            substrLength = hhsearchConfig.parameters['substrLength']):

        self.query = query
        self.inputDirectory = inputDirectory
        self.user = user
        self.hhsearchPATH = hhsearchPATH
        self.hhsuiteDir = hhsuiteDir
        self.uniprotDir = uniprotDir
        self.hhsearchDB = hhsearchDB
        self.hhblits = hhsearchConfig.parameters['hhblits']     #hhsuiteDir+'/bin/hhblits'
        self.hhmake  = hhsearchConfig.parameters['hhmake']      #hhsuiteDir+'/bin/hhmake'
        self.hhsearch = hhsearchConfig.parameters['hhsearch']       #hhsuiteDir+'/bin/hhsearch'
        self.addss  = hhsearchConfig.parameters['addss.pl']     #Path to run addss.pl
        self.runMake = runMake
        self.runSearch = runSearch
        self.runComplex = runComplex
        self.runMonomer = runMonomer
        self.runGenome  = runGenome
        self.seqfileExtension = seqfileExtension
        self.tmpDir = hhsearchConfig.parameters['tmpDir']+ '/'+user+'/hhsearch_'+os.path.basename(query)+'/'#Defined in HHsearchMain
        self.substrPosition = substrPosition
        self.substrLength = substrLength

    def HHsearchMain(self):
        """
        Main function for HHsearch Wrappers.  Creates temporary directories, and
        runs hhsearch and hhmake
        """
        commandfunctions.mkdirs(self.tmpDir)
        os.chdir(self.tmpDir)
        queryList = self.GetQueryList()
        for query in queryList:
            fastaPath,hhmPath,hhrPath,workDir = self.GetAndMakePaths(query)
            self.RunHHmake_HHsearch(query,fastaPath,hhmPath,hhrPath,workDir)
        commandfunctions.RemoveDirectory(self.tmpDir)

    def GetQueryList(self):
        """
        Creates queryList either from a newline delimited file or string inputed
        from command line. 

        Returns:
            queryList ([str]): List of query strings ie 117eA or 117eA-117eB.
        Notes:
            The directory structure is different for a single chain suchas 
            117eA compared to 117eA-117eB.  Read readme for futher details.
        """
        queryList = []
        if os.path.isfile(self.query):
            f = open(self.query)
            queryList = [line.strip() for line in f]
            f.close()
        else:
            queryList = [self.query.strip()]
        return queryList


    def RunHHmake_HHsearch(self,query,fastaPath,hhmPath,hhrPath,workdir):
        """
        Creates working directories for HHsearch and HHmake.  Will
        run either/both HHmake and HHsearch depending on if self.runMake
        and self.runSearch are set to true.

        Arguments:
            query (str): query string ie 117eA or 117eA-117eB
            fastaPath (str): Full Path to location of fasta file.
            hhmPath (str). Full Path to location of hhm file 
            hhrPath (str). Full Path to location of hhr file
            workDir (str). Path to root temporary directory.
        Notes:
            If the fasta file is /home/bgovi/test/117eA.fasta and
            the hhm file is /home/bgovi/test/117eA.hhm and the 
            hhr file is /home/bgovi/test/117eA.hhr,
            fastaPath == hhmPath == hhrPath == /home/bgovi/test/

        """ 
        if '-' in query and self.runComplex:
            queryA,queryB = query.split('-')
            workhhmA = workdir + '/hhmA/' 
            workhhmB = workdir + '/hhmB/'
            workhhrA = workdir + '/hhrA/'
            workhhrB = workdir + '/hhrB/'
            if self.runMake:
                commandfunctions.mkdir(workhhmA)
                commandfunctions.mkdir(workhhmB)
                self.Make(queryA,fastaPath,hhmPath,workhhmA)
                self.Make(queryB,fastaPath,hhmPath,workhhmB)
            if self.runSearch:
                commandfunctions.mkdir(workhhrA)
                commandfunctions.mkdir(workhhrB)
                self.Search(queryA,hhmPath,hhrPath,workhhrA)
                self.Search(queryB,hhmPath,hhrPath,workhhrB)
        else:
            workhhm = workdir +'/hhm/'
            workhhr = workdir +'/hhr/'
            if self.runMake:
                commandfunctions.mkdir(workhhm)
                self.Make(query,fastaPath,hhmPath,workhhm)
            if self.runSearch:
                commandfunctions.mkdir(workhhr)
                self.Search(query,hhmPath,hhrPath,workhhr)

    def GetAndMakePaths(self, query):
        """
        Creates input and ouptput paths for running hhsearch/hhmake for
        single chains and complexes.  A complex is two query chains concatenated
        with a '-' such as 117eA-117eB.  The complex and single chain have 
        seperate directory structures.  Refer to readme for more information.

        Returns:
                        query (str): query string ie 117eA or 117eA-117eB
                        fastaPath (str): Full Path to location of fasta file.
                        hhmPath (str). Full Path to location of hhm file 
                        hhrPath (str). Full Path to location of hhr file
                        workDir (str). Path to root temporary directory.
                Notes:
                        If the fasta file is /home/bgovi/test/117eA.fasta and
                        the hhm file is /home/bgovi/test/117eA.hhm and the 
                        hhr file is /home/bgovi/test/117eA.hhr,
                        fastaPath == hhmPath == hhrPath == /home/bgovi/test/
        """
        workDir = self.tmpDir+query+'/'
        commandfunctions.mkdir(workDir)
        if '-' in query and self.runComplex:
            fastaPath = self.inputDirectory+'/'+query +'/'
            hhmPath = fastaPath
            hhrPath = fastaPath
            commandfunctions.mkdir(fastaPath)
            return fastaPath, hhmPath, hhrPath,workDir

        elif self.runGenome:
            startSubstr = int(self.substrPosition)
            endSubstr = startSubstr + int(self.substrLength)
            commandfunctions.mkdir(self.inputDirectory+'/hhm/')
            commandfunctions.mkdir(self.inputDirectory+'/hhr/')
            fastaPath = self.inputDirectory+'/fasta/'+query[startSubstr:endSubstr]+'/'
            hhmPath = self.inputDirectory+'/hhm/'+query[startSubstr:endSubstr]+'/'
            hhrPath = self.inputDirectory+'/hhr/'+query[startSubstr:endSubstr]+'/'
            commandfunctions.mkdir(fastaPath)
            commandfunctions.mkdir(hhmPath)
            commandfunctions.mkdir(hhrPath)
            return fastaPath, hhmPath, hhrPath,workDir
        elif self.runMonomer:
            fastaPath = self.inputDirectory+'/'+query +'/'
            hhmPath = fastaPath
            hhrPath = fastaPath
            return fastaPath, hhmPath, hhrPath,workDir

        else:
            if self.runComplex:
                errorString = "runComplex set to true but query name is in wrong format. "
                errorString+="Make sure query has name chain1-chain2 ie 12asA-12asB."
                raise ValueError(errorString)
            else:
                errorString = "runComplex, runMonomer and runGenome all set to false.  Make sure "
                errorString += " one of the parameters is set to true"
                raise ValueError(errorString)
            
    def Make(self,queryName, fastaPath, hhmPath,workingDir):
        """
        Make runs several programs from the hhsuite library inorder to 
        generate a profile hidden markov model (hhm) 
        Arguments
            queryName (str): Name of target sequence. For example 117eA, which
                is chain A from the pdbfile 117e
            fastaPath (str): Path to location of query sequence file.  If
                sequence file is located in /home/bgovi/117eA.seq. 
                queryPath = /home/bgovi/
            hhmPath (str): Path to output location of hhm file.
            workingDir (str): Temporary directory used for running hhsuite 
                functions. ie. /tmp/bgovi/117eA/ 
        ConfigArguments:
            self.seqFileExtenstion (str,optional): Sequence Files may have 
                different extensions or no extension.  If 117eA is sequence file
                then seqfileExtension == '', if 117eA.seq is sequence file then 
                seqfileExtension == '.seq'
        Output
            Makes an hhm file for the query sequence, and outputs it to the 
            queryPath. For example /queryPath/queryName.hhm
        """
        #if results file exists, skip
        if os.path.isfile(hhmPath+'/'+queryName+".hhm"):
            print(queryName,"hhm already completed")
            return

        #copy sequence file to workingDir
        fastaFile = fastaPath+'/'+queryName+self.seqfileExtension
        queryFilePath = workingDir +'/'+queryName + self.seqfileExtension
        shutil.copyfile(fastaFile, queryFilePath)

        #run hhblits
        oa3mFile = workingDir+'/'+queryName +".a3m"
        cmd = [self.hhblits, '-i', queryFilePath, '-oa3m', oa3mFile]
        cmd = self._AppendArguments(cmd, hhsearchConfig.hhblitsParameters)
        stdout, stderr = commandfunctions.Execute(cmd,raiseStdError = False)
        #run addss.pl
        cmd = [self.addss, oa3mFile]
        stdout, stderr = commandfunctions.Execute(cmd,raiseStdError=False,
                            environment = os.environ)
        #run hhmake
        cmd = [self.hhmake,'-i',oa3mFile]
        cmd = self._AppendArguments(cmd, hhsearchConfig.hhmakeParameters)
        stdout, stderr = commandfunctions.Execute(cmd,raiseStdError=False)
        shutil.copyfile(workingDir+'/'+queryName+'.hhm',
                        hhmPath+'/'+queryName+'.hhm')

    def Search(self,queryName,hhmPath,hhrPath, workingDir):
        """
        Runs HHsearch and creates hhr results file.

        Requires:
            queryName (str): Name ot target hhm.  For example 117eA which is 
                chainA from the pdbfile 117e.  The hhm file = 117eA.hhm
            hhmPath (str): Full path to location of hhm file.
            hhrPath (str): Full path to output hhr file.
            workingDir (str): Temporary directory to run HHsearch
        Outputs:
            Copies the queryName.hhr results file to hhrPath
        """
        #check if hhr file completed
        if os.path.isfile(hhrPath+'/'+queryName+'.hhr'):
            print(queryName, 'hhr already completed')
            return  
        #copy hhm file to workingDir
        shutil.copyfile(hhmPath+'/'+queryName+'.hhm',
                        workingDir+'/'+queryName+'.hhm')
        #run hhsearch
        cmd = [self.hhsearch,'-i',workingDir+'/'+queryName+'.hhm']
        cmd = self._AppendArguments(cmd, hhsearchConfig.hhsearchParameters)
        stdout,stderr = commandfunctions.Execute(cmd, raiseStdError=False)
        shutil.copyfile(workingDir+'/'+queryName+'.hhr',
                        hhrPath+'/'+queryName+'.hhr')

    def Submit(self):
        """
        For Job Submission.  Job parameters can be changed in hhsearchConfig.py
        """
        walltime = hhsearchConfig.jobParameters['walltime']
        priority = hhsearchConfig.jobParameters['priority']
        maxJobs =  hhsearchConfig.jobParameters['maxJobs']
        forceSubmit = hhsearchConfig.jobParameters['forceSubmit']
        hhsearchJob = JOBS(self.user,walltime=walltime,priority=priority,
                        maxJobs=maxJobs)
        jobBaseName = 'hhsearch_'
        scriptName = os.path.basename(__file__)
        if self.runMake and not self.runSearch:
            jobBaseName = 'hhmake_'


        if not os.path.isfile(self.query):
            jobInput = [hhsearchConfig.pythonPath,
                       hhsearchConfig.parameters['hhsearchPATH']+'/'+scriptName,
                       '-q',self.query,'-iDir',self.inputDirectory]
            self._AppendArguments_toJobInput(jobInput)
            outputDir = self._JobOutputDir(self.query)
            jobOutput = self._JobOutputList([self.query])
            jobName = jobBaseName+os.path.basename(self.query)
            hhsearchJob.SubmitJob(jobInput, jobOutput, jobName, outputDir,
                                    ForceSubmit=forceSubmit)
        else:
            queryList = self.GetQueryList()
            if (hhsearchConfig.jobParameters['batchSize'] <= 1 or
                len(queryList) == 1):
                for query in queryList:
                    jobInput = [hhsearchConfig.pythonPath,
                       hhsearchConfig.parameters['hhsearchPATH']+'/'+scriptName,
                       '-q',query,'-iDir',self.inputDirectory]
                    self._AppendArguments_toJobInput(jobInput)
                    if self.runGenome:
                        outputDir = self.inputDirectory+'/record/'
                    else:
                        outputDir = self._JobOutputDir(query)
                    jobOutput = self._JobOutputList([query])
                    jobName = jobBaseName+query
                    hhsearchJob.SubmitJob(jobInput, jobOutput, jobName,
                                            outputDir,ForceSubmit=forceSubmit)
            else: #write files
                i = 0
                batchsize = hhsearchConfig.jobParameters['batchSize']
                if not isinstance(batchsize, int):
                    raise ValueError("batchSize must be an integer in jobParameters in hhsearchConfig.py")
                for j in xrange(0,len(queryList),batchsize):
                    tmpqueryList = queryList[j:j+batchsize]
                    outputDir = self.inputDirectory+'/record/'
                    commandfunctions.mkdir(outputDir)
                    flistPath = self._WriteQueryList(i,tmpqueryList,outputDir)
                    jobInput = [hhsearchConfig.pythonPath,
                       hhsearchConfig.parameters['hhsearchPATH']+'/'+scriptName,
                       '-q',flistPath,'-iDir',self.inputDirectory]
                    self._AppendArguments_toJobInput(jobInput)
                    jobOutput = self._JobOutputList(tmpqueryList)
                    jobName = jobBaseName+os.path.basename(flistPath)
                    hhsearchJob.SubmitJob(jobInput, jobOutput,jobName,outputDir,
                                            ForceSubmit=forceSubmit)
                    i+=1

    def _AppendArguments_toJobInput(self,jobInput):
        """
        Checks to see if runMake,runComplex or runSearch are true.  If true 
        append to job submission
        """
        if self.runMake:
            jobInput.append('-make')
        if self.runSearch:
            jobInput.append('-search')
        if self.runComplex:
            jobInput.append('-complex')
        if self.runGenome:
            jobInput.append('-genome')
        if self.runMonomer:
            jobInput.append('-monomer')

                
    def _WriteQueryList(self,jobNum, queryList,outputDir):
        """
        Write Query List for batch submission and return name and full path of 
        list file. For file with a long list of query submissions.  This 
        recieves a subset of a list and write it to a file for job submission.

        Arguments:
            jobNum (int): Non negative integer for job submission.  Used for 
                providing unique jobNames during batch submission.  ie 
                submitting groups of query targets as one job.
            queryList ([str]): List of jobs query targets for submission.
            outputDir ([str]): Directory where list of query targets will be 
                written to
        Outputs:
            list_jobNum (file): Newline delimieted list of query targets.  if 
                    jobNum = 1
                    the output file is list_1
        Returns:
            foutName (str): Full path name of output list file. 
        """
        foutName = outputDir+'/list_'+str(jobNum)
        fout = open(foutName,'w')
        for query in queryList:
            fout.write(query+'\n')
        fout.close()
        return foutName
            


    def _JobOutputDir(self,query):
        """
        Returns directory where job submission scripts are located for qsub 
        command.
        """
        if '-' in self.query:
            outputDir  = self.inputDirectory +'/'+query+'/record/'
        else:
            outputDir = self.inputDirectory+'/record/'
        return outputDir

    def _JobOutputList(self,queryList):
        """
        Returns a long list of files that are generated if a job has already 
        been submitted and completed. Prevents unintentional reruns of already 
        completed searches.
        """
        outputList = []
        for query in queryList:
            if '-' not in query:
                startPos = int(hhsearchConfig.parameters['substrPosition'])
                endPos = startPos + int(hhsearchConfig.parameters['substrLength'])
                if self.runMake == True:
                    hhmPath = self.inputDirectory+'/hhm/'+ query[startPos:endPos]+'/'+query+'.hhm'
                    outputList.append(hhmPath)
                if self.runSearch == True:
                    hhrPath = self.inputDirectory+'/hhr/'+ query[startPos:endPos]+'/'+query+'.hhr'
                    outputList.append(hhrPath)
            else:
                queryAB = query.split('-')
                if self.runMake == True:
                    hhmPathA = self.inputDirectory+'/'+ query+'/'+queryAB[0]+'.hhm'
                    hhmPathB = self.inputDirectory+'/'+ query+'/'+queryAB[1]+'.hhm'
                    outputList.append(hhmPathA)
                    outputList.append(hhmPathB)
                if self.runSearch == True:
                    hhrPathA = self.inputDirectory+'/'+ query+'/'+queryAB[0]+'.hhr'
                    hhrPathB = self.inputDirectory+'/'+ query+'/'+queryAB[1]+'.hhr'
        return outputList


    def _AppendArguments(self,cmd, argumentDictionary):
        """
        Appends flag,alguments in argumentDictionary to command list (cmd).
        The dictionaries are from the hhsearchConfig.py file
        Arguments:
            cmd ([str]): Starting command list for running system process
            arguementDictionary ({}): Dictionary of additional commands. where the
                key is the command argument flag and the value is the argument for
                each flag. 
        Returns:
            cmd ([str]): Appends key,value pairs from argument dictionary to the original
                    cmd list.
        """
        for flag, argument in argumentDictionary.iteritems():
            if argument != '':
                cmd.append(flag)
                cmd.append(argument)
        return cmd  

def ParseCommandLine():
    """
    Parses command line for arguments for running hhsearch.py

    Returns:
        args (Namespace from argparse): This namespace will contain all the 
        parameters from the command line
    """
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__),
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    arguments = parser.add_argument_group(title="Arguments")
    queryHelp = "query file or directory ie -q 117eA-117eB or -q "\
            "/home/username/querylist.txt where querylist.txt is a newline "\
            "delimited list of targets ie 117eA-117eB 12asA-12asB etc. or "\
            "117eA 117eB 12asA etc. The full file path is required if "\
            "inputing a list."
    arguments.add_argument('-q',action="store",required=True,help=queryHelp)

    inputDirHelp = "The directory that contains the target directory data. "\
        "For complexes like 117eA-117eB if the directory 117eA-117eB full "\
        "path is /home/user/targets/117eA-117eB the inputDirectory should "\
        "store /home/user/targets/ ie -iDir /home/user/targets/ For single "\
        "chain targets. the input directory should contain a fasta, hhm and "\
        "hhr directory.  if these directories are in /home/user/targets/  "\
        "the -iDir = /home/user/targets/"
    arguments.add_argument('-iDir',action="store",required=True,
        help=inputDirHelp)

    jobHelp="Job flag is used for job submission, default runs program locally."
    arguments.add_argument('-job',action="store_true",default=False,
        help=jobHelp)

    makeHelp = "Optional flag for running hhmake.  Defualt operation is "\
        "determined by runMake in hhsearchConfig.py, which can be set to "\
        "true and this flag can be ignored."
    arguments.add_argument('-make',action="store_true",
        default=hhsearchConfig.parameters['runMake'],help=makeHelp)

    searchHelp = "Optional flag for running hhsearch.  Default operation is "\
        "determined by runSearch in hhsearcConfig.py, which can be set to "\
        "true and this flag can be ignored."
    arguments.add_argument('-search',action='store_true',
        default=hhsearchConfig.parameters['runSearch'],help=searchHelp)

    complexHelp ="optional flag required for running complexes. ie targets "\
        "in form 117eA-117eB.  The hyphen seperates the two chains.  "\
        "Default operation is determined by runComplex in hhsearchConfig.py"
    arguments.add_argument('-complex',action='store_true',
        default=hhsearchConfig.parameters['runComplex'],help=complexHelp)

    monomerHelp = "optional flag required for running monomers. ie targets "\
        "in form 117eA.  The hyphen seperates the two chains.  Default "\
        "operation is determined by runComplex in hhsearchConfig.py"
    arguments.add_argument('-monomer',action='store_true',
        default=hhsearchConfig.parameters['runMonomer'],help=monomerHelp)

    genomeHelp = "optional flag required for running monomers in the genome "\
        "wide format.  The hyphen seperates the two chains.  Default "\
        "operation is determined by runComplex in hhsearchConfig.py"
 
    arguments.add_argument('-genome',action='store_true',
        default=hhsearchConfig.parameters['runGenome'],help=genomeHelp)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = ParseCommandLine()
    query = args.q
    inputDirectory = args.iDir
    submitJob = args.job
    runMake = args.make
    runSearch = args.search
    runComplex = args.complex
    runMonomer = args.monomer
    runGenome  = args.genome

    countTrue = 0
    for i in [runComplex,runMonomer,runGenome]:
        if i:
            countTrue+=1

    if countTrue == 0:
        print("Flag error.  runComplex, runMonomer and runGenome are all ",
                "set to False")
        print("One of the flags ie -complex,-monomer or -genome are needed")
        print("Or runComplex, runMonomer or runGenome can be set to True in ",
                "hhsearchConfig.py")
        sys.exit() 

    if countTrue > 1:
        print("More than one of the flags ie -complex,-monomer,-genome are ",
                "used or ")
        print("the default settings of runComplex,runMonomer and runGenome ",
                "have more")
        print("than one set to true. Only one of the formats are allowed per ",
                "run.")
        print("Check hhsearchConfig.py or make sure only one of the flags is ",
                "used at the command line")
        sys.exit()

    hhsearch = HHSEARCH(query, inputDirectory,runMake=runMake,
            runSearch=runSearch,runComplex = runComplex,runMonomer=runMonomer,
            runGenome=runGenome)
    if hhsearch.runMake == False and hhsearch.runSearch == False:
        print("Both runSearcn and runMake are set to false, causing this ",
                "script to due nothing.")
        print("Either set runSearch or runMake to True in hhsearchConfig.py ",  
                "or use the -make or -search flags")
        sys.exit() 

    if submitJob:   #submit job
        hhsearch.Submit()
    else:   #run locally
        hhsearch.HHsearchMain() 
