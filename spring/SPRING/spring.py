#!/usr/bin/env python
"""
The spring module for threading of protein complexes.  

All default command line inputs and paths can be changed in
springConfig.py 
"""
#python 2.x and 3.x cross compatibility
from __future__ import print_function, division
import sys
if sys.version_info[0] >= 3:
    xrange = range

#python libraries
import os
import numpy as np
from operator import itemgetter
import sys
import shutil
import argparse

#user libraries
import springConfig
ROOT_PATH = springConfig.parameters['rootDir']
sys.path.append( springConfig.parameters['pdbPath'])
from pdbchains import PDBchains, ReadFasta
from alignment import TMalign, TMscore, NWalign
from interface import InterfaceContacts, FractionNativeContacts
sys.path.append( springConfig.parameters['jobPath'])
from jobsubmit import JOBS
sys.path.append( springConfig.parameters['hhsearchPath'])
from hhfunctions import HHSEARCH_Functions
from hhsearch import HHSEARCH
sys.path.append( springConfig.parameters['commandPath'])
from commandfunctions import Execute, mkdirs, mkdir, RemoveDirectory

class SPRING:
    """
    This class contains the data structures and functions for running SPRING 
    complex threading.

    Attributes:
        RefCoreToChain: A dictionary that chains in monomer library to 
                        consituents in the complex library
        BioMol: A dicitionary containing all the names of complex interactions 
                in the PDB
        query: The target ie 12asA-12asB or the full path of a file containing a
               list of targets.
        indexFile: The path to the index file
        dfireFile: The path to the dfire file:
        DfireEnergy: The stored dfire energy obtained from dfireFile
        W0: First Weight for spring score
        W1: Second Weight for spring score
        Normalize: A value to adjust the spring score to make values over 10 
                   is significant.
        seqidCutoff: Sequence id threshold hold for excluding templates
        user: username ie 'bgovi'
        AminoToCode: Converts amino character to interger list.  Integer list 
                     used for dfire energy
        inputDirectory: path to location of target directories
        outputDirectory: Opitonal. default is the target directory
        run_hhsearch: Runs hhsearch if hhr files not present in target directory
        workDir: Temporary directroy created to run spring.
        springOutputDirectory: outputDirectory+'/SPRING/'. Results are sent here
        queryDirectory: Full path to target directory where all input files are 
                        located.
        queryName: Name of target ie 12asA-12asB
        queryA: Name of first chain ie 12asA
        queryB: Name of second chain ie 12asB
        sequenceA: Sequence of first chain
        sequenceB: Sequence of second chain
        fileModelChainA: Optionaly path for monomer model or template used for 
                         superimposing onto complex framework for chainA
        fileModelChainB: Optionaly path for monomer model or template used for 
                         superimposing onto complex framework for chainB
        hhrFileA: Path to the hhr results file for chain A ie 12asA.hhr
        hhrFileB: Path to the hhr results file for chain B ie 12asB.hhr
        qtemplatesA: Stored information in hhr file for chainA into pandas data 
                     structure
        qtemplatesB: Stored information in hhr file for chainB into pandas data 
                     structure
        monomerA: PDBchains object for top ranked hhr result or chain from 
                  fileModelChainA 
        monomerB: PDBchains object for top ranked hhr result or chain from 
                  fileModelChainB 
        dockMonomer: If true SPING outputs models with top ranked monomer
                    superimposed onto dimer framework

    Methods:
        RunSpring: Main function for running spring.
        SwitchTopMonomer: Uses pdb from fileModelChain if present else use pdb 
                          from hhr file
        SpringTemplateSearch: Main function for searching for complex templates
        Submit: Submits spring jobs
        _ExtendJobInput: Adds arguments to spring jobs
        RemoveMonomerFilePaths: Clears temporary data paths for next iteration 
                                of spring.
        UpdateQueryData:Updates temporary data paths for nex iteration of spring
        GetQueryPath: Checks for paths containing hhr files and fileModelChainA.
        GetQueryList: Gets list of targets from query if its a file.
        SearchIndex: Finds initial list of interactions from monomer search and 
                     index file 
        MakeModelAndScore: Makes and stores Template Models and Scores
        MakeModel: Makes template comlex models.
        ChainOrderSwitch: If queryA and queryB are not homodimers. Searches 
                          database by swapping A and B
        ScoreTemplate: Calculate Dfire
        StoreDfire: Stores dfireFile into DfireEnergy
        DFIRE: Calculates DFIRE potential for complex model
        HHsearchResultsComplete: Checks if hhsearch results complete. Can run 
                                 hhsearch if paraemter in config file selected.
        ReadHHsearchResults: Reads hhearch results file and stores information 
                             into qtemplatesA and qtemplatesB
        StoreIndex: Stores index file into RefCoreToChain and BioMol
        WriteOutput: Write out init.dat TemplateSummary.txt and Templates to 
                     springOutputDirectory
        Completed: Check if spring has been ran and completed.
    """
    def __init__(self,query,inputDirectory,user=springConfig.parameters['user'],
        outputDirectory='',W0=springConfig.parameters['W0'], 
        W1=springConfig.parameters['W1'],
        Normalize=springConfig.parameters['Normalize'],
        seqidCutoff = springConfig.parameters['seqidCutoff'], 
        fileModelChainA='',fileModelChainB='',hhrFileA='',hhrFileB='',
        runHHsearch = springConfig.parameters['runHHsearch'],
        indexFile = springConfig.parameters['indexFile'], 
        dfireFile = springConfig.parameters['dfireFile'],
        dockMonomer = False):

        self.RefCoreToChain = {}
        self.BioMol = {}
        self.query = query
        self.indexFile = indexFile
        self.dfireFile = dfireFile
        self.DfireEnergy = []
        self.W0 = W0
        self.W1 = W1
        self.Normalize = Normalize
        self.seqidCutoff = seqidCutoff
        self.user = user
        #for dfire. everything stored as one big array
        self.AminoToCode={'A':0,'C':1,'D':2,'E':3,'F':4,'G':5,'H':6,'I':7,'K':8,
           'L':9,'M':10,'N':11,'P':12,'Q':13,'R':14,'S':15,'T':16,'V':17,'W':18,
           'Y':19,'B':20,'Z':20,'X':20}
        self.query = query
        self.inputDirectory = inputDirectory
        self.outputDirectory = outputDirectory
        self.run_hhsearch = runHHsearch
        self.workDir = springConfig.parameters['tmpDir']+'/'+ \
                        self.user+'/SPRING_'+os.path.basename(self.query)+'/'
        self.dockMonomer = dockMonomer
        #These attributes will change for each queryName if query is a list of targets
        self.springOutputDirectory = ''
        self.queryDirectory = ''
        self.queryName = ''
        self.queryA = ''
        self.queryB = ''
        self.sequenceA = ''
        self.sequenceB = ''
        self.fileModelChainA = fileModelChainA
        self.fileModelChainB = fileModelChainB
        self.hhrFileA = hhrFileA
        self.hhrFileB = hhrFileB
        self.qtemplatesA = None
        self.qtemplatesB = None
        self.monomerA = None
        self.monomerB = None

    def __repr__(self):
        return "SPRING("+self.query+")"
    def __str__(self):
        return "member of SPRING class"

    def RunSpring(self):
        """
        Main function for running spring.py
        1. Create and change to temporary working directory
        2. Store index.txt file
        3. Store dfire.txt file
        4. Itereate spring search through all queries in the queryList
            1. Clear Paths from Previous Run
            2. Update Paths for new query
            3. Read In query sequences
            4. Read and/or run hhsearch for both query chains
            5. Check if alternate top monomer present for spring docking
            6. Search for complex templates.
            7. Write output TemplateSummary.txt, init.dat and Templates to
                spring output directory
        """
        mkdirs(self.workDir)
        os.chdir(self.workDir) 
        print("Store Index")
        self.StoreIndex(self.indexFile)
        print("Store Dfire")
        self.StoreDfire(self.dfireFile)
        PDBfunctions = PDBchains()
        queryList = self.GetQueryList(self.query)
        for query in queryList:
            self.RemoveMonomerFilePaths(len(queryList)) #clears user input if running spring in batch mode
            print(query)
            #store queryName, queryA, queryB, hhr files, and monomerA, monomerB locations
            self.UpdateQueryData(query)
            if self.Completed(self.springOutputDirectory): #check if arleady completed
                continue
            print("Read Fasta")
            self.ReadInSequences(self.queryA,self.queryB,self.queryDirectory)
            #check if hhsearch complete. If not complete, if run_hhearch true its ran else skip query
            if not self.HHsearchResultsComplete():
                continue
            print("Read HHsearch")
            print("Using HHsearch Files ",self.hhrFileA, self.hhrFileB)
            self.ReadHHsearchResults(self.hhrFileA, self.hhrFileB)
            templateDir = self.springOutputDirectory+'/Templates/'
            modelDirectory = self.springOutputDirectory+'/Models/'
            initDir = self.springOutputDirectory
            if self.qtemplatesA.empty or self.qtemplatesB.empty:
                self.WriteOutput([],modelDirectory,self.queryName,templateDir,
                        initDir) #send file with DONE
                continue
            self.SwitchTopMonomer() #if user provides monomer path switch top monomer
            templateList = self.SpringTemplateSearch(self.qtemplatesA,
                                self.qtemplatesB)
            print("Writting Templates and Template Summary to File")
            self.WriteOutput(templateList,modelDirectory,
                    self.queryName,templateDir,initDir)
        RemoveDirectory(self.workDir)


    def SwitchTopMonomer(self):
        """
        If user provided top monomers present in springOutputDirectory. Use them
        instead of top rank monomer from hhsearch for docking.
        
        Notes:
            The user provide monomer has a specific format for example the 
            target 117eA-117eB.  Inside the springOutputDirectory there can be 
            two files 117eA_template.pdb and 117eB_template.pdb or the user can 
            provide a unique path at instantiation of object or on the command 
            line.  The command line option is ignored if running in batch mode. 
        Modifies:
            self.monomerA (pdbchains): Top ranked monomer from hhsearch or user 
                        provided monomer for chain A
            self.monomerB (pdbchains): Top ranked monomer from hhsearch or user 
                        provided monomer for chain B
        """
        chainA = self.fileModelChainA
        if os.path.exists( chainA ):
            self.monomerA = PDBchains()
            self.monomerA.Read(chainA, ca_only=True, allow_alt_loc = False,
                                allow_insert = False)
        chainB = self.fileModelChainB
        if os.path.exists( chainB ):
            self.monomerB = PDBchains()
            self.monomerB.Read(chainB, ca_only=True, allow_alt_loc = False, 
                                allow_insert = False)
            

    def SpringTemplateSearch(self,qtemplatesA,qtemplatesB,
            homodimerThreshold = springConfig.parameters['homodimerThreshold']):
        """
        Main function for spring complex threading.

        qtemplatesA and qtemplatesB: (pandas dataframe): Data structure obtained
            from /source/hhsearch/hhfunctions.py qtemplates stores alignment and
            summary data for hhr files for both chains of query Template
        Returns:
            templateList ([ springscore,[coreName,pdbCore], [partnerName,pdbBP],
                energyScores]):  Template ranked/sorted by springscore. coreName
                and partnerName are the names of both chains for the complex 
                template.  pdbCore and pdbBP are the models for the complex 
                template stored using pdbchains.  energyscores are the various 
                sequence alignment scores from hhsearch. 
        """
        templateList = []
        for i in xrange(0,2):
            if i ==1: 
                nwout = NWalign(self.sequenceA, self.sequenceB)
                if (nwout['seqidA'] >= homodimerThreshold and
                        nwout['seqidB'] >= homodimerThreshold):
                    continue
                else: #switch search order for heterodimer. Now seqB is core and seqA is partner
                    qtemp = qtemplatesA
                    qtemplatesA = qtemplatesB
                    qtemplatesB = qtemp 
            print("Search Index")
            complexList = self.SearchIndex(qtemplatesA,qtemplatesB)
            tmpNum = 0
            for a in xrange(0,len(complexList)):
                tmpNum += len(complexList[a][5])
            print("Number of Possible Templates:", tmpNum)
            print("Make Template Models and obtain Spring Score")

            if i == 0:
                templateList.extend( self.MakeModelAndScore(complexList,i) )
            else: #switch rank1s back
                tmpList2 = []
                tmpList2.extend(self.MakeModelAndScore(complexList,i) )
                templateListTmp = []
                for j in xrange(0,len(tmpList2)): #switch core and partner order
                    #templateList ([ springscore, [core,pdbCore], 
                    #   [partnerName,pdbBP],energyScores] )
                    #energy scores ([dfireEnergy,tmscoreA,tmscoreB,coreZscore,
                    #   coreProbability,coreEvalue,coreSeqid, 
                    #   partnerZscore,partnerProbability,partnerEvalue,
                    #   partnerSeqid])   
                    springscore = tmpList2[j][0]
                    core = tmpList2[j][1]
                    partner = tmpList2[j][2]
                    energyList = tmpList2[j][3]
                    partnerScores = energyList[7:]
                    coreScores = energyList[3:7]
                    dfire = energyList[0]
                    tmscoreA = energyList[1]
                    tmscoreB = energyList[2]
                    energy = [dfire,tmscoreB,tmscoreA]
                    energy.extend(partnerScores)
                    energy.extend(coreScores)
                    templateListTmp.append([springscore,partner, core,energy])
                templateList.extend(templateListTmp)
        print("Rank Templates")
        if not templateList: #list empty
            return templateList
        templateList = sorted(templateList, key=itemgetter(0), reverse=True)
        return templateList

    def Submit(self,recordDir = '', 
        walltime = springConfig.jobParameters['walltime'], 
        priority = springConfig.jobParameters['priority'],
        maxJobs =  springConfig.jobParameters['maxJobs'],
        forceSubmit = springConfig.jobParameters['forceSubmit']):
        """
        Function for submiting spring script
        """
        springJob = JOBS(self.user,walltime=walltime,priority=priority,
                        maxJobs=maxJobs)
        #scriptName = os.path.basename(__file__)
        scriptName = springConfig.parameters['springPY']
        jobInput = [springConfig.pythonPath,
                   springConfig.parameters['springPath']+'/'+scriptName,
                   '-q',self.query,'-iDir',self.inputDirectory]
        self._ExtendJobInput(jobInput)
        if recordDir == '':
            outputDir = self.inputDirectory+'/record/'
            mkdir(outputDir)
        else:
            outputDir = recordDir
            mkdir(outputDir)
        jobName = "SPRING_"+os.path.basename(self.query)
        jobOutput = []
        jobOutput.append(self.springOutputDirectory+'/TemplateSummary.txt')
        springJob.SubmitJob(jobInput,jobOutput,jobName,outputDir,
                            ForceSubmit=forceSubmit)

    def JobFailed(self,target,recordDir=''):
        """
        Checks to see if spring job failed. Returns False if job completed.
        Returns true if pbs error and output files exists but spring output
        files dont
        """
        jobName = "SPRING_"+os.path.basename(self.query)
        springJob = JOBS(self.user)
        outputDir = recordDir
        if recordDir == '':
            if not os.path.isfile(self.query):
                outputDir = self.inputDirectory+'/record/'
                mkdir(outputDir)
            else: 
                outputDir = self.inputDirectory+'/'+self.query+'/'
                mkdir(outputDir)
                outputDir+='/record/'
                mkdir(outputDir)
        if self.outputDirectory != '':
            self.springOutputDirectory = self.outputDirectory+'/SPRING/'
        else:
            self.springOutputDirectory = self.inputDirectory+'/'+target+'/SPRING/'
        jobOutputList = [self.springOutputDirectory+'/TemplateSummary.txt',self.springOutputDirectory+'/init.dat']
        if self.dockMonomer:
            modelPath = self.springOutputDirectory+'/Models/model1.pdb'
            jobOutputList.append(modelPath)
        return springJob.JobFailed(outputDir,jobName,jobOutputList)

    def _ExtendJobInput(self,jobInput):
        """
        Add arguments to job script if provided by user
        """
        if self.outputDirectory != '':
            jobInput.extend(['-oDir',self.outputDirectory])
        if self.fileModelChainA != '':
            jobInput.extend(['-c1',self.fileModelChainA])
        if self.fileModelChainB != '':
            jobInput.extend(['-c2',self.fileModelChainB])
        if self.hhrFileA != '':
            jobInput.extend(['-hhr1',self.hhrFileA])
        if self.hhrFileB != '':
            jobInput.extend(['-hhr2',self.hhrFileB])
        if self.run_hhsearch:
            jobInput.append('-hhsearch')
        if self.W0 != springConfig.parameters['W0']:
            jobInput.extend(['-W0','%.2f'%self.W0])
        if self.W1 != springConfig.parameters['W1']:
            jobInput.extend(['-W1','%.2f'%self.W1])
        if self.Normalize != springConfig.parameters['Normalize']:
            jobInput.extend(['-Norm','%.2f'%self.Normalize])
        if self.dockMonomer:
            jobInput.append('-dockMono')

    def RemoveMonomerFilePaths(self,lenQueryList):
        """
        If query list larger than 1 only default paths for hhr files and
        model templates are used.  User input paths for these files are
        ignored.
        Argument:
            lenQueryList (int): Number of query in queryList.  queryList
                can be greater than 1 if self.query is a file of target
                complexes
        """
        if lenQueryList > 1:
            self.fileModelChainA = ''
            self.fileModelChainB = ''
            self.hhrFileA =''
            self.hhrFileB =''

    def UpdateQueryData(self,query):
        """
        Updates paths regarding target query.

        Arguments:
            query (str): query target ie 117eA-117eB
        Modifies:
            self.queryName (str): Name of query target ie 117eA-117eB
            self.queryA    (str): Name of 1st chain of query ie 117eA
            self.queryB    (str): Name of 2nd chain of query ie 117eB
            self.queryDir  (str): Full path of query output directory 
            self.hhrFileA  (str): Full path to hhr file for chain A
            self.hhrFileB  (str): Full path to hhr file for chain B
            self.fileModelChainA (str): Full path to template model for chain A
            self.fileModelChainB (str): Full path to template model for chain B
        """
        #store queryName, queryA,queryB and Sequences
        self.queryA = query.split('-')[0]
        self.queryB = query.split('-')[1]
        self.queryName = query
        self.queryDirectory = self.inputDirectory+'/'+query+'/'
        if not os.path.exists(self.queryDirectory):
            mkdirs(self.queryDirectory)
        if self.outputDirectory != '':
            mkdirs(outputDirectory)
            self.springOutputDirectory = self.outputDirectory+'/SPRING/'
        else:
            self.springOutputDirectory = self.queryDirectory+'/SPRING/'
        mkdirs(self.springOutputDirectory)

        if self.fileModelChainA == '':
            self.fileModelChainA = self.GetQueryPath(self.springOutputDirectory,
                                self.queryDirectory,self.queryA,'_template.pdb')
        if self.fileModelChainB == '':
            self.fileModelChainB = self.GetQueryPath(self.springOutputDirectory,
                                self.queryDirectory,self.queryB,'_template.pdb')
        if self.hhrFileA == '':
            self.hhrFileA = self.GetQueryPath(self.springOutputDirectory,
                                self.queryDirectory,self.queryA,'.hhr')
        if self.hhrFileB == '':
            self.hhrFileB = self.GetQueryPath(self.springOutputDirectory,
                                self.queryDirectory,self.queryB,'.hhr')
    def GetQueryPath(self,firstPath,backupPath,queryName,extension):
        """
        Checks for file path in springOutputDirectory and inputDirectory.
        Returns file path in query inputDirectory if it exists and doesn't exist
        in springOutptuDirectory.  Else it returns default file path to 
        springOutputDirectory
        """
        path1 = firstPath +queryName+extension
        path2 = backupPath+queryName+extension
        path3 = path1
        path4 = path1
        if extension == '.hhr':
            i = springConfig.PDBparameters['substrPosition']
            j = springConfig.PDBparameters['substrEnd']
            path3 = self.inputDirectory+'/hhr/'+queryName[i:j]+'/'
            path3 += queryName+extension
            path4 = self.inputDirectory+'/../hhr/'+queryName[i:j]+'/'
            path4 += queryName+extension
        if os.path.exists(path2) and not os.path.exists(path1):
            return path2
        elif os.path.exists(path1):
            return path1
        else:
            if os.path.exists(path3):
                return path3
            if os.path.exists(path4):
                return path4
            else:
                return path1


    def GetQueryList(self,query):
        """
        Generates query list from query if query is a newline delimited file
        containing a list of targets:

        Arguments:
            query (str or file path): query is either a target ie 117eA-117eB or
                the full path to a file containing list of targets.
        Returns:
            queryList ([str]): List of targets. ie 
                ['117eA-117eB','12asA-12asB',...]

        """
        if os.path.isfile(self.query):
            f = open(self.query)
            queryList = [line.strip() for line in f]
            f.close()
            return queryList
        else:
            queryList = [self.query]
            return queryList

    def SearchIndex(self,qtemplatesA, qtemplatesB):
        """
        Returns a list of interactions that are found using the reference
        lookup search.

        Arguments:
            qtemplatesA and qtemplatesB (pandas data frame): The HHsearch
                results for queryA and queryB. Obtained from
                HHSEARCH_Functions.ReadHHsearchSummary. This moduled is found
                in the source/hhsearch directory inside hhfunctions.py .  
                The variables have sequence to template alignments and 
                similarity scores.
        Returns:
            complexList =[[core,coreZscore,coreProbability,coreEvalue,coreSeqid,
                        [partner,partnerZscore,partnerProbability,
                        partnerEvalue, partnerSeqId]] ]
            complexList are all complex templates found in index who 
                have core reference chains which are in qtemplatesA and 
                binding partners references in qtemplatesB.  
                The data contained in core and partner are described
                in the StoreIndex function that is a part of the 
                spring module.
        Notes:
            The similarity scores contained in complexList are from the 
            hhr results file from running hhsearch
        """
        Btemplates = {}
        coreList = []
        complexList = [] #[ [core,partner,coreZscore,partnerZscore] ]
        for i in xrange(0,qtemplatesB.shape[0]):
            reference = qtemplatesB.iloc[i]['Template']
            refBzscore = qtemplatesB.iloc[i]['SpringZscore']
            refBProbability = qtemplatesB.iloc[i]['Probability']
            refBE_value = qtemplatesB.iloc[i]['Evalue']
            refBseqid   = qtemplatesB.iloc[i]['SeqID']
            values = (refBzscore, refBProbability,refBE_value,refBseqid)
            Btemplates[reference] = values

        for i in xrange(0,qtemplatesA.shape[0]):
            reference = qtemplatesA.iloc[i]['Template']
            refAzscore = qtemplatesA.iloc[i]['SpringZscore']
            refAProbability = qtemplatesA.iloc[i]['Probability']
            refAE_value = qtemplatesA.iloc[i]['Evalue']
            refAseqid   = qtemplatesA.iloc[i]['SeqID']
            values = [reference, refAzscore, refAProbability,refAE_value,
                        refAseqid]
            if reference in self.RefCoreToChain:
                for coreChain in self.RefCoreToChain[reference]:
                    coreList.append(values)
        if not coreList:
            return complexList

        coreHit = {}
        complexListCutoff = []
        for coreData in coreList:
            reference = coreData[0]
            refAzscore = coreData[1]
            refAProbability = coreData[2]
            refAE_value = coreData[3]
            refAseqid = coreData[4]
            for coreChain in self.RefCoreToChain[reference]:
                if coreChain in coreHit:
                    continue
                print(coreChain)
                coreHit[coreChain] = None
                if coreChain not in self.BioMol:
                    continue
                for biomol in self.BioMol[coreChain].keys():
                    interactions = self.BioMol[coreChain][biomol]
                    core = interactions[0]
                    partners = interactions[1]
                    numberInteractions = interactions[2]
                    partnerList = []
                    partnerListCutoff = []
                    for partner in partners:
                        partnerCore = partner[0]
                        partnerBioMol_CoreOrPartner_chainNum = partner[1]
                        partnerRef = partner[2]
                        if partnerRef in Btemplates:
                            values = Btemplates[partnerRef]
                            refBzscore, refBProbability,refBE_value,refBseqid = values
                            partnerList.append( [partner, refBzscore,
                                    refBProbability,refBE_value,refBseqid] )
                            if refBE_value <= 0.01:
                                partnerListCutoff.append( [partner, refBzscore,
                                    refBProbability,refBE_value,refBseqid] )
                    complexList.append( [core,refAzscore,refAProbability,
                                            refAE_value,refAseqid,partnerList] )
                    if refAE_value <= 0.01 and len(partnerListCutoff) > 0:
                        complexListCutoff.append( [core,refAzscore,refAProbability,
                            refAE_value,refAseqid,partnerListCutoff])

        if len(complexListCutoff) > 1:
            return complexListCutoff
        else:
            return complexList

    def MakeModelAndScore(self,complexList, iteration):
        """
        Make and store monomers

        Arguments:
            complexList = [ [coreName,coreZscore,coreProbability,coreEvalue,
                        coreSeqid,[partnerName,partnerZscore,partnerProbability,
                        partnerEvalue, partnerSeqId]] ]
                complexList are all complex templates found in index who 
                have core reference chains which are in qtemplatesA and 
                binding partners references in qtemplatesB.  
                The data contained in core and partner are described
                in the StoreIndex function that is a part of the 
                spring module.
            iteration (int) 0 or 1. 0 for first order 1 for second.
        Returns:
            templateList ([ springscore, [core,pdbCore], [partnerName,pdbBP],
                            energyScores] ):
                springscore (float): Overall spring score
                core (str): Name of core complex template pdb
                pdbCore (pdbchains): Template model for the first chain
                partnerName (str): Name of binding partner to core pdb
                pdbBP (pdbchains): Template model for the second chain
                energy scores ([dfireEnergy,tmscoreA,tmscoreB,coreZscore,
                            coreProbability,coreEvalue,coreSeqid, 
                            partnerZscore,partnerProbability,partnerEvalue,
                            partnerSeqid])   
        """
        templateList = []
        monomerA,chainALen,monomerB,chainBLen = self.ChainOrderSwitch(iteration)
        count = 0
        for interaction in complexList:
            coreName = interaction[0]
            coreZscore = interaction[1]
            coreProbability = interaction[2]
            coreEvalue = interaction[3]
            coreSeqid = interaction[4]
            partners = interaction[5]
            if len(partners) == 0:
                continue
            chainA_name = coreName[0]+'/'+coreName[1]
            pdbCore, tmscoreA = self.MakeModel(chainA_name, monomerA, chainALen)
            for partner in partners:
                print(count)
                partnerName = partner[0]
                partnerZscore = partner[1]
                partnerProbability = partner[2]
                partnerEvalue = partner[3]
                partnerSeqid = partner[4]
                chainB_name = partnerName[0]+'/'+partnerName[1]
                pdbBP, tmscoreB = self.MakeModel(chainB_name,monomerB,chainBLen)
                dfireEnergy = self.DFIRE(pdbCore, pdbBP)
                springscore = self.ScoreTemplate(coreZscore,partnerZscore,
                                    tmscoreA,tmscoreB,dfireEnergy)
                energyScores = [dfireEnergy,tmscoreA,tmscoreB,coreZscore,
                                coreProbability,coreEvalue,coreSeqid,
                                partnerZscore, partnerProbability, 
                                partnerEvalue, partnerSeqid]
                templateList.append([springscore, [coreName,pdbCore], 
                                        [partner[0],pdbBP],energyScores] )
                count += 1

        return templateList

    def MakeModel(self,templateName, monomer, chainLen,
        start = springConfig.PDBparameters['substrPosition'],
        end = springConfig.PDBparameters['substrEnd'],
        pdbPath = springConfig.PDBparameters['pdbDir'] + '/',
        pdbExt = springConfig.PDBparameters['pdbFileExtension'] ):
        """
        Map monomer sequence alignment to complex template framework by TMalign.

        Arguments:
            templateName (str): Name of template in index.txt ie 1o1mB/0_0_0
            monomer (pdbchains): Top ranked monomer from hhsearch or user
                provided monomer
            chainLen (int): Lenght of query chain 1 or chain 2
            start   (int): Start of slice for template Name. 
                        (This is part of PDB structure read README for info)
            end     (int): End of slice for template Name. 
            pdbPath (str): Path to complex PDB directory
            pdbExt  (str): Extension of pdbfile normally '.pdb'
        Returns:
            modelChain (pdbchains): Monomer alignment mapped onto complex 
                                    template chain
            tmscore    (float):  TMscore from TMalign renormalized by chainLen
        """
        templateFile =  pdbPath+templateName[start:end]+'/'+templateName+pdbExt
        cplxChain = PDBchains()
        cplxChain.Read(templateFile, ca_only=True, allow_alt_loc = False, 
                        allow_insert = False)
        transfromMatrix, tmout = TMalign(monomer, cplxChain)
        monomerAlignment = tmout['seqModel_align']
        cplxChainAlign  = tmout['seqNative_align']
        tmscore =  (tmout['TMscore_normNative']*len(cplxChain.sequence[0]) )/(1.0*chainLen)
        cplxResNum  = []
        cplxResName = []
        monomerResNum = []
        monomerResName = []
        cplxSeq = []
        monoPos = 0
        cplxPos = 0
        for i in xrange(0,len(monomerAlignment) ):
            if monomerAlignment[i] != '-' and cplxChainAlign[i] != '-':
                monoCApos = monomer.ca_pos[0][monoPos]
                cplxCApos = cplxChain.ca_pos[0][cplxPos]
                monomerResNum.append(  monomer.atom_info[0][monoCApos].res_num )
                monomerResName.append( monomer.atom_info[0][monoCApos].res_name)
                cplxResNum.append(  cplxChain.atom_info[0][cplxCApos].res_num )
                cplxResName.append( cplxChain.atom_info[0][cplxCApos].res_name )
                monoPos += 1
                cplxPos += 1
            elif cplxChainAlign[i] != '-':
                cplxPos += 1
            else:
                monoPos += 1
        modelChain = cplxChain.SliceResNum(0,cplxResNum, ca_only = True)
        #replace sequence and residue numbers
        for i in xrange(0,len(modelChain.sequence[0] ) ):
            modelChain.atom_info[0][i].res_num  = monomerResNum[i]
            modelChain.atom_info[0][i].res_name = monomerResName[i]
            modelChain.atom_info[0][i].temp_num = cplxResNum[i]
            modelChain.atom_info[0][i].temp_res = cplxResName[i]
        return modelChain, tmscore  
    
    def ChainOrderSwitch(self, iteration):
        """
        Returns order for monomer and sequences.  Spring searches library using
        both chains in both orders.
        """
        if iteration == 0 or iteration == 1:
            pass
        else:
            raise ValueError(' iteration must be 0 or 1\n')
        if iteration == 0:
            chainALen = len(self.sequenceA)*1.0
            chainBLen = len(self.sequenceB)*1.0
            monomerA = self.monomerA
            monomerB = self.monomerB
            return monomerA, chainALen, monomerB, chainBLen
        else:#for switching orders
            chainALen = len(self.sequenceB)*1.0
            chainBLen = len(self.sequenceA)*1.0
            monomerA = self.monomerB
            monomerB = self.monomerA
            return monomerA, chainALen, monomerB, chainBLen



    def ScoreTemplate(self,zscoreA,zscoreB,tmscoreA,tmscoreB,dfireEnergy):
        """
        Calculates Spring Score
        """

        minZscore = min(zscoreA,zscoreB)
        minTMscore = min(tmscoreA,tmscoreB)
        springscore = minZscore + self.W0*minTMscore + self.W1*dfireEnergy
        springscore = springscore/self.Normalize
        return springscore

    def StoreDfire(self,DfireFile):
        """
        Stores dfire.txt into DfireEnergy array

        The Dfire Energy array is the dfire pairwise distance potential 
        for all 20 standard amino acids.  The structure of the array is 
        as follows.

        PairWisePotential=ResToCode['Amino']*21*20 + AminoToCode['Amino']*20 + 
        distBin  
    
        ResToCode is a dicitionary that converts one letter Amino Acid codes 
        into integers.  This is done in alphabitical order so A = 0, C=1,D=2 and
        so on. Aminos designated with B,Z, or X are given a value of 20.
        
        The 21 represents the number of types which is 21.  The standard 20
        plus one for other.  The 20 is the number of bins used in the dfire
        potential.  The 20 bins range from 0-10 angstroms with 0.5 Angstrom
        intervals.

        distBin is the bined pairwise distance.  distBin = int(2*distance)
            anything greater than 20 should not be sent to this array
        
        Arguments:
            DfireFile (str): Full file path of the dfire potential  
        """
        f = open(DfireFile)
        for eval in f:
            eval = float(eval.strip())
            self.DfireEnergy.append(eval)
        f.close()
    
    def DFIRE(self,pdbCore, pdbBP,chainA=0,chainB=0):
        """
        Calculates the DFIRE Potential. 

        For the structure of the dfire potential read the StoreDfire 
        docstring/comments.

        Arguments:
            pdbCore and pdbBP (pdbchains object): Each object contains one chain
                of the protein dimer.  
        Returns:
            dfireEnergy (float): Pairwise dfire energy potential of dimer 
                                 interface.
        """
        dfireScale = 2.0 #for bin size 0-10 with 0.5 step size bins
        distMax = 10
        dfireEnergy = 0
        intcontactOut = InterfaceContacts(pdbCore,pdbBP, chainA=chainA,
                                chainB=chainB,distanceCutoff = distMax)
        contactMatrix = intcontactOut[0]
        distanceMatrix = intcontactOut[1]
        #get array slice of row if true.
        distBins = distanceMatrix*dfireScale
        distBins.astype(int)
        #I,J are the row and column indexes where there is an interface contact.
        #So the first contact is at I[0],J[0]
        I,J   = np.where(contactMatrix)
        #calculate dfire energy
        for x in xrange(0,len(I)):
            i = I[x]
            j = J[x]            
            aminoBP = pdbBP.sequence[chainA][j]
            aminoCore = pdbCore.sequence[chainB][i] 
            distBin = distBins[i][j]
            try:
                index = self.AminoToCode[aminoBP]*21*20 + self.AminoToCode[aminoCore]*20 + distBin
                dfireEnergy += self.DfireEnergy[int(index)]
            except:
                pass
        return dfireEnergy


    def DFIRE_Monomer(self,pdbCore,chainNum=0):
        """
        Calculates the DFIRE Potential. 

        For the structure of the dfire potential read the StoreDfire 
        docstring/comments.

        Arguments:
            pdbCore (pdbchains object): Protein Chain  
        Returns:
            dfireEnergy (float): Pairwise dfire energy potential of protein fold
        """
        dfireScale = 2.0 #for bin size 0-10 with 0.5 step size bins
        distMax = 10
        dfireEnergy = 0
        intcontactOut = InterfaceContacts(pdbCore,pdbCore, chainA=chainNum,
                                chainB=chainNum,distanceCutoff = distMax)
        contactMatrix = intcontactOut[0]
        distanceMatrix = intcontactOut[1]
        #get array slice of row if true.
        distBins = distanceMatrix*dfireScale
        distBins.astype(int)
        #I,J are the row and column indexes where there is an interface contact.
        #So the first contact is at I[0],J[0]
        I,J   = np.where(contactMatrix)
        #calculate dfire energy
        for x in xrange(0,len(I)):
            i = I[x]
            j = J[x]
            if i > j:
                continue
            aminoBP = pdbCore.sequence[chainNum][j]
            aminoCore = pdbCore.sequence[chainNum][i] 
            distBin = distBins[i][j]
            try:
                index = self.AminoToCode[aminoBP]*21*20 + self.AminoToCode[aminoCore]*20 + distBin
                dfireEnergy += self.DfireEnergy[int(index)]
            except:
                pass
        return dfireEnergy

    def HHsearchResultsComplete(self):
        """
        Check existence of hhr results files for both chains of the query.
        If they are absent and run_hhsearch set to true hhsearch is ran on the 
        query chain
        
        Returns:
            True if both hhr files present False otherwise
        """
        hhrFiles = [self.hhrFileA, self.hhrFileB]
        outStr = ''
        i = 0
        for hhrFile in hhrFiles:
            if not os.path.exists(hhrFile):
                outStr+='For query '+self.queryName+"hhr file: "+hhrFile+" missing\n"
                if self.run_hhsearch:
                    fastaPath = self.queryDirectory
                    hhmPath   = self.springOutputDirectory
                    hhrPath   = self.springOutputDirectory
                    query = ''
                    if i == 0:
                        query = self.queryA
                    else:
                        query = self.queryB
                    hhsearch = HHSEARCH(query,self.user,runMake=True,
                            runSearch=True,runMonomer=True, 
                            seqfileExtension=springConfig.parameters[
                                                    'seqFileExtension'])
                    hhsearch.RunHHmake_HHsearch(query,fastaPath,hhmPath,hhrPath,
                            self.workDir)
                i+=1
                    
        if outStr != '' and not self.run_hhsearch:
            print(outStr,end='')
            print("HHR files missing for query "+self.queryName,
                " and run_hhsearch set to False, either change",
                " runHHmake to True at")
            print("the command line or springConfig.py or run HHsearch",
                 " seperately then rerun SPRING")
            return False

        if (not os.path.exists(self.hhrFileA) or 
            not os.path.exists(self.hhrFileB)):
            print("HHR files missing for query "+self.queryName,
                 " and HHsearch failed to generate hhr files.",
                 "  Please check to see")
            print("why hhsearch failed then rerun SPRING")
            return False
        
        else:
            return True

    def ReadHHsearchResults(self,hhsearchResultsA, hhsearchResultsB,
                pdbFileExt = springConfig.PDBparameters['monomerPDBExtension'],
                writeMonomer = True):
        """
        Reads the HHsearch file for each query sequence and stores the results.
        The results are filtered based on the sequence identity cutoff.

        Arguments:

        Modifies:
            self.qtemplatesA and self.qtemplatesB (pandas data frame): Contains
                template names and alignment scores for each 
                template.  For further detais of whats contained
                check ReadHHsearch in pdbchains module.
            self.monomerA and self.monomerB (pdbchains): template models from
                1st rank hhsearch results model
        Output:
            Writes 1st rank template models to springOutputDirectory as 
            rank1_chain1.pdb and rank1_chain2.pdb
        """
        hhr1 = HHSEARCH_Functions(hhsearchResultsA, GetAlignmentData = True)
        hhr2 = HHSEARCH_Functions(hhsearchResultsB, GetAlignmentData = True)
        resultsA = hhr1.summary
        resultsB = hhr2.summary

        if resultsA.empty or resultsB.empty:
            self.qtemplatesA = resultsA
            self.qtemplatesB = resultsB
            return
        #sequence identity cutoff
        if self.seqidCutoff >= 1 or self.seqidCutoff < 0:
            self.qtemplatesA = resultsA
            self.qtemplatesB = resultsB
        else:
            self.qtemplatesA = resultsA[ resultsA['SeqID'] <= self.seqidCutoff ]
            self.qtemplatesB = resultsB[ resultsB['SeqID'] <= self.seqidCutoff ]
        if self.qtemplatesA.empty or self.qtemplatesB.empty:
                return
        rank1A = self.qtemplatesA.iloc[0]['Template']
        rank1B = self.qtemplatesB.iloc[0]['Template']
        monomerPDB = springConfig.PDBparameters['monomerpdbDir']
        start = springConfig.PDBparameters['substrPosition']
        end =   springConfig.PDBparameters['substrLength'] + start
        templateA_Path = monomerPDB+'/'+rank1A[start:end]+'/'+rank1A+pdbFileExt
        templateB_Path = monomerPDB+'/'+rank1B[start:end]+'/'+rank1B+pdbFileExt
        self.monomerA = hhr1.MakeModel(templateA_Path, rank1A)
        self.monomerB = hhr2.MakeModel(templateB_Path, rank1B)
        #renumber monomerB
        chain1Len = len(self.sequenceA)
        for i in xrange(0,len(self.monomerB.atom_info[0]) ):
            self.monomerB.atom_info[0][i].res_num += chain1Len

        if writeMonomer:
            self.monomerA.Write(self.springOutputDirectory+'/rank1_chain1.pdb',
                            write_template = True)
            self.monomerB.Write(self.springOutputDirectory+'/rank1_chain2.pdb',
                            write_template = True)

 
    def StoreIndex(self,indexFile):
        """
        This function stores the spring index file.

        File Format for index.txt
        
        Example segment from index.txt file             
        column1  column2   column3
        12e8H    0_0_0     12e8L
        12e8H    0_1_2     12e8L
        12e8H    0_1_3     12e8H
        137lA    0_0_0     137lA
        137lA    0_1_1     212lA
        137lA    1_0_0     137lA
        137lA    1_1_1     212lA
        12e8L    1_0_0     12e8L
        12e8L    1_1_1     12e8H

        The index.txt has three columns:

        Column 1 is the core pdb chain.  Every row that 
        has the same core pdb chain refers to an interaction
        involving that chain.

        Column2: contains three pieces of information.
        biomolecule_CoreOrPartner_chainNumber
            
        biomolecule refers to which biomolecule the interaction
        comes from.  Each pdbfile has one or several biomoleulces. 
        Each biomolecule is a different set of possilby biologically 
        relevant orientations of each macromolecule with respect to 
        one another.  For example one PDB file may contain a dimer
        with two biomolecules.  Each biomolecule representing a different
        interface.  As you can from the example above 11gsA only has
        one biomolecule designated with a 0 where as 137lA has two 
        biomolecules designated with a 0 and 1.

        CoreOrPartner refers to a row being a core or a partner.  For
        example 12e8H has three rows.  The first row with 0_0_0 was 
        the chain selected to check for interactions (the core).  
        Every macromolecule in contact with it is considered to be
        a partner.  Here the core is designated with a 0.  This
        core has two binding partners.  Each binding partner is
        designated with the number 1.

        chainNumber refers to the order of chains in a pdb file.
        The first chain is 0 the second is 1 and so on.  In the
        12e8H rows the pdb file contains 4 chains, but only the
        last two are involved in interactions with the first chain.
        So there is no interaction of the core with the second
        chain in the PDB file.

        Column 3 contains the reference pdb.  Each pair of sequences
        is searched as monomers through a monomeric database containing
        pdbs found in the reference column. For example if a query
        sequence matches to 12e8L all of the rows containing 12e8H
        and 12e8L in the first column would be used as for further 
        evaluation as a complex template because the both have core
        pdbs that share the 12e8L as a reference.  Each reference
        in column 3 shares a high degree of sequence and structural
        similarity to the pdb chain in column 1.   

        Modifies:  
            RefCoreToChain (dictionary): For each reference a list
                of all pdbchains that it references are stored.
                The list is actually another dictionary to allow
                for faster comparision. 
                Example Structure:
                RefCoreToChain[referencePDB] = {core1:None,core2:None,...}
            bioMolecule (dictionary): For each core a list (stored as a
                dictionary), of each biomolecules is stored. Each biomolecule 
                points to a list containing three values [core,[partners],
                number interactions] The core is the core pdb stored as a tuple
                containing the row information found in index.txt.  The second 
                value is a list of all partners.  Each partner is stored as a 
                tuple containing its respective row data from index.txt.  
                Finally the third value is number of interactions for the core.
                Example Structure:
                bioMolecule[core][biomoleculeNumber] = [chain,[partner],
                                                        numInteractions]
                bioMolecule['12e8L'][1] = [(12e8L,1_0_0,12e8L),
                                            [(12e8L,1_1_1,12e8H)],1]
        """
        f = open(indexFile)
        for line in f:
            tmp = line.split()
            coreChain = tmp[0]
            chainInfo = tmp[1] #biomolecule_core/partner_chainNum
            repChain = tmp[2]
            chain = (coreChain,chainInfo,repChain)

            tmp = chainInfo.split('_')
            bioMolecule = int(tmp[0])
            CoreOrPartner = int(tmp[1])
            chainNum = int(tmp[2])
            #store reference -> to all pdb chain it references  
            if repChain not in self.RefCoreToChain and CoreOrPartner == 0:
                self.RefCoreToChain[repChain] = {}

            if CoreOrPartner == 0:
                if coreChain not in self.RefCoreToChain[repChain]:
                    self.RefCoreToChain[repChain][coreChain] = None


            if coreChain not in self.BioMol:
                if CoreOrPartner == 0:
                    self.BioMol[coreChain] = {bioMolecule:[chain,[],0]}
                else:
                    self.BioMol[coreChain] = {bioMolecule:['',[chain],1]}   
            else:
                if bioMolecule in self.BioMol[coreChain]:
                    if CoreOrPartner == 0:
                        self.BioMol[coreChain][bioMolecule][0] = chain
                    else:
                        self.BioMol[coreChain][bioMolecule][1].append(chain)
                        self.BioMol[coreChain][bioMolecule][2] += 1
                else:
                    if CoreOrPartner == 0:
                        self.BioMol[coreChain][bioMolecule] = [chain,[],0]
                    else:
                        self.BioMol[coreChain][bioMolecule] = ['',[chain],1]
        f.close()

    def WriteOutput(self,templateList,modelDirectory,query,templateDir, initDir,
            maxNumTempl =springConfig.parameters['maxNumberTemplates'],
            maxNumInit = springConfig.parameters['maxNumberTemplateInit'],
            pdbPath = springConfig.PDBparameters['pdbDir'],
            subStart = springConfig.PDBparameters['substrPosition'],
            subEnd = springConfig.PDBparameters['substrLength'],
            pdbExt = springConfig.PDBparameters['pdbFileExtension'],
            maxModels = springConfig.parameters['maxModels']):
        """
        Writes output TemplateSummary.txt init.dat (for TACOS simulation)  and 
        Templates to File.

        Arguments:
            templateList ([ springscore, [core,pdbCore], [partnerName,pdbBP],
                            energyScores] ):
                springscore (float): Overall spring score
                core (str): Name of core complex template pdb
                pdbCore (pdbchains): Template model for the first chain
                partnerName (str): Name of binding partner to core pdb
                pdbBP (pdbchains): Template model for the second chain
                energy scores ([dfireEnergy,tmscoreA,tmscoreB,coreZscore,
                   coreProbability,coreEvalue,coreSeqid, 
                   partnerZscore,partnerProbability,partnerEvalue,partnerSeqid])
            query (str): Name of query ie 117eA-117eB
            templateDir (str): Directory where templates models are written
            initDir     (str): Directory where init.dat is written
            manNumTempl (int): Maximum number of templates to be written
            maxNumInit  (int): Maximum number of templates in init.dat
        """
        modelComplex = PDBchains()
        if self.dockMonomer:
            self.ReadHHsearchResults(self.hhrFileA,self.hhrFileB,writeMonomer=False)
            self.SwitchTopMonomer()
            self.monomerA.ChangeChainID(0,'A')
            self.monomerB.ChangeChainID(0,'B')
            modelComplex.Append(self.monomerA,chains=[0])
            modelComplex.Append(self.monomerB,chains=[0])
            modelComplex.Renumber(renumber_res=False,start_reset=False)
            mkdir(modelDirectory)

        mkdir(initDir)
        mkdir(templateDir)
        complexLen = len(self.sequenceA) + len(self.sequenceB)
        numTemplates = 0
        fout = open(initDir+'/TemplateSummary.txt','w')
        fout.write(' '.join(['%-5s'%'#','%-25s'%'Complex Template',
            '%12s'%'SpringScore','%9s'%'LenAligned','%6s'%'Dfire',
            '%9s'%'Coverage','%6s'%'SeqidA','%6s'%'TMscoreA', '%6s'%'ZscoreA',
            '%6s'%'ProbabilityA', '%6s'%'EvauleA','%6s'%'SeqidB',
            '%6s'%'TMscoreB', '%6s'%'ZscoreB','%6s'%'ProbabilityB',
            '%6s'%'EvauleB']))
        fout.write('\n')
        init = PDBchains()
        numberTemplates = 0
        writeSummaryAll = springConfig.parameters['WriteFullSummary']
        for i in xrange(0,len(templateList) ):
            tmpNum = i+1
            compleX = PDBchains()
            templateData = templateList[i]
            springscore  = templateData[0]
            coreName, corePDB = templateData[1]
            coreName_Full = coreName[0]+'/'+coreName[1]
            partnerName, partnerPDB = templateData[2]
            partnerName_Full = partnerName[0]+'/'+partnerName[1]
            dfire,tmA,tmB,zsA,probA,evalA,seqidA,zsB,probB,evalB,seqidB = templateData[3]
            compleX.Append(corePDB, chains = [0])
            compleX.Append(partnerPDB, chains = [0])
            if tmpNum <= maxNumTempl:
                compleX.Write(templateDir+'/Template'+str(tmpNum)+'.pdb')
            LenAligned = len(compleX.sequence[0])+len(compleX.sequence[1])
            coverage = (1.0*LenAligned)/complexLen
            complexTemplateName = coreName_Full+'-'+partnerName_Full

            corePath = pdbPath+'/'+coreName_Full[subStart:subStart+subEnd]+'/'+coreName_Full
            partnerPath = pdbPath+'/'+partnerName_Full[subStart:subStart+subEnd]+'/'+partnerName_Full

            #need full paths to corePDB and templatePDB
            seqidA = self.TemplateQuerySeqID(corePath,self.sequenceA,
                                            pdbFileExt = pdbExt)
            seqidB = self.TemplateQuerySeqID(partnerPath,self.sequenceB,
                                            pdbFileExt = pdbExt)
            seqidAB = (seqidA+seqidB)/2.0
            #need global seqid.  Need to run nwalign.  Need
            fout.write(' '.join(['%-5s'%tmpNum,'%-25s'%complexTemplateName,
                    '%12.3f'%springscore,'%10s'%LenAligned,'%6.3f'%dfire,
                    '%9.3f'%coverage,'%6.3f'%seqidA,'%8.3f'%tmA, '%7.3f'%zsA,
                    '%6.3f'%probA, '%11.3g'%evalA,'%7.3f'%seqidB,'%6.3f'%tmB,
                    '%9.3f'%zsB,'%7.3f'%probB, '%11.3g'%evalB]))
            fout.write('\n')
            if tmpNum <= maxNumInit:
                init.Append(compleX,chains=[0])
                init.Extend(i,compleX,1,useRenumber = False)
                zscoreINIT = 5.0*min(tmA,tmB)
                init_tempHeader = ''.join(['%5s'%LenAligned,'%8.3f'%zscoreINIT,
                            '%4s'%tmpNum,'%27s'%complexTemplateName,
                            '%8.3f'%seqidAB,
                            '%8.3f'%coverage+"(="+str(LenAligned)+'/'+str(complexLen)+')',
                            " (L_ali,Z,i,pdb,id,cov)"])
                init.template_headers[i] = init_tempHeader
                numberTemplates += 1

            if self.dockMonomer and tmpNum <= maxModels:
                tmOut = TMalign(modelComplex,compleX,modelchain=0,nativechain=0,
                            applyTransform=True) 
                tmOut = TMalign(modelComplex,compleX,modelchain=1,nativechain=1,
                            applyTransform=True)
                modelOut = modelDirectory+'/model'+str(tmpNum)+'.pdb'
                modelComplex.Write(modelOut,write_connect=True)
            if not writeSummaryAll:
                if tmpNum > maxNumTempl and tmpNum > maxNumInit:
                    break
        fout.write('DONE')
        fout.close()
        if init.num_chains > 0:
            initFile = initDir +'/init.dat'
            init.init_header = ''.join(['%5s'%numberTemplates,'%5s'%complexLen,
                             '%5s'%len(self.sequenceA), " (N_temp, Lch, Lch1)"])
            init.Write(initFile,write_init=True)
             

    def TemplateQuerySeqID(self,template,querySeq,
                pdbFileExt = springConfig.PDBparameters['pdbFileExtension']):
        """
        Aligns original template sequence to query sequence and returns sequence
        identity
        Arguments:
            template (str): Template Name ie 3p8tA/0_0_0
        """
        tmpPDB = PDBchains()
        tmpPDB.Read(template+pdbFileExt,ca_only=True)
        tmpSeq = tmpPDB.sequence[0]
        nwout = NWalign(tmpSeq,querySeq)
        seqid = nwout['seqidB']
        return seqid

    def ReadInSequences(self, queryA, queryB, queryDirectory):
        """
        Reads in query sequences.
        Arguments:
            queryA   (str): Name of query for chain A
            queryB   (str): Name of query for chain B
            queryDir (str): Path to query sequences
        Modifies:
            self.sequenceA and self.sequenceB by storing query sequence
        """
        seqFileExt = springConfig.parameters['seqFileExtension']
        seqAFile = queryDirectory+queryA+seqFileExt
        seqBFile = queryDirectory+queryB+seqFileExt
        fastaPath1 = self.inputDirectory+'/fasta/'
        fastaPath2 = self.inputDirectory+'/../fasta/'
        start = springConfig.PDBparameters['substrPosition']
        end   = springConfig.PDBparameters['substrEnd']
        seqExt= springConfig.PDBparameters['seqFileExtension']
        if not os.path.exists(seqAFile):
            print('Warning',seqAFile,'does not exist try alternate file paths')
            seqAFile = fastaPath1+queryA[start:end]+'/'+queryA+seqExt
            if not os.path.exists(seqAFile):    
                seqAFile = fastaPath2+queryA[start:end]+'/'+queryA+seqExt
            print('Using new sequence path',seqAFile)
        if not os.path.exists(seqBFile):
            print('Warning',seqBFile,'does not exist try alternate file path')
            seqBFile = fastaPath1+queryB[start:end]+'/'+queryB+seqExt
            if not os.path.exists(seqBFile):    
                seqBFile = fastaPath2+queryB[start:end]+'/'+queryB+seqExt
            print('Using new sequence path',seqBFile)
        headers, sequences = ReadFasta(seqAFile)
        self.sequenceA = sequences[0]
        headers, sequences = ReadFasta(seqBFile)
        self.sequenceB = sequences[0]

    def Completed(self, queryDir='',verbose=True):
        """
        Checks if summary.txt has been completed (Looks for line starting
        with DONE)
    
        Arguments:
            queryDir (str): Full path to query directory
            verbose  (bool): If true prints if target completed
        Returns:
            bool: True if DONE in summary.txt, False otherwise
        """
        if queryDir =='':
            queryDir = self.springOutputDirectory
        fileName = queryDir + "TemplateSummary.txt"
        if not os.path.exists(fileName):
            return False

        f = open(fileName)
        count = 0
        for line in f:
            if line[0:4] == "DONE":
                f.close()
                out = queryDir + " search completed"
                if verbose:
                    print(out)
                return True
            count+=1
        f.close()

        if self.dockMonomer:
            modelPath = queryDir+'/Models/model1.pdb'
            if not os.path.exists(modelPath) and count < 2:
                return False
        return False


    def _pdbFilePath(self,pdbDirectory, indexRow):
        """
        Returns full path name of pdb template in splits directory
        index Row comes from index.txt
        12e8H    0_0_0     12e8L

        indexRow = (12e8H,0_0_0,12e8L)
        """
        
        pdbFile = pdbDirectory +'/'+ indexRow[0][0:2] + '/' +indexRow[0] +'/'+ indexRow[1]
        pdbFile += ".pdb"
        return pdbFile

    def _rankFilePath(self, chainDirectory, chainName):
        """
        Returns full path name of pdb template contained in the 
        HHsearch database.
        """
        f = chainDirectory + '/' +chainName[0:2] + '/' + chainName
        return f
#############################################################################
# The functions below are meant to be used by the main function or as an    #
# interface for another module that wishes to run spring as a function while#
# using the command line interface.  Other methods for running spring would #
# be to import the class or run as a subprocess                             #
#############################################################################

def ParseCommandLine():
    """
    Parses command line for arguments for running SPRING.py
    
    Returns:
        args (Namespace for argparse): This namespace will contain all the 
        parameters from the command line
    """
    defaultParam = springConfig.parameters
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__),
                description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    #Required Arguments
    arguments = parser.add_argument_group(title='Arguments')

    queryHelp = "query file or directory ie -q 117eA-117eB or -q "\
            "/home/username/querylist.txt where querylist.txt is a newline"\
            " delimited list of targets ie 117eA-117eB 12asA-12asB etc. "\
            "The full file path is required if inputing a list.  If the query "\
            "is a list the -c1 or -c2 flags are ignored.  HHsearch is ran "\
            "unless the the default names for the monomer templates are "\
            "present."
    arguments.add_argument('-q',action="store",required=True,help=queryHelp)

    inputDirHelp = "The directory that contains the target directory data. "\
            "For complexes like 117eA-117eB if the directory 117eA-117eB full"\
            "  path is /home/user/targets/117eA-117eB the inputDirectory "\
            "should store /home/user/targets/ ie -iDir /home/user/targets/."
    arguments.add_argument('-iDir',action="store",required=True,
            help=inputDirHelp)

    jobHelp = "Job flag is used for job submission, default runs program "\
            "locally."
    arguments.add_argument('-job',action="store_true",default=False,
              help=jobHelp)

    outputDirHelp = "Optional. Full path for COTH output.  Default is the "\
            "targetDirectory which is the inputDirectory/query/ is "\
            "/home/username/117eA-117eB/"
    arguments.add_argument('-oDir',action="store",default='',help=inputDirHelp)

    chainAfileHelp ="Optional.  Full path to chain 1 template obtained from "\
            "monomer threading."
    arguments.add_argument('-c1',action="store",default='',help=chainAfileHelp)

    chainBfileHelp="Optional. Full path to chain 2 template obtained from "\
            "monomer threading."
    arguments.add_argument('-c2',action="store",default='',help=chainBfileHelp)

    hhrFile1Help = "Optional. Full path to the hhr file results for chain A."\
            "Default checks for existence in inputDir/target/SPRING/ then "\
            "inputDir/target/ ie /home/user/117eA-117eB/SPRING/117eA.hhr or "\
            "/home/user/117eA-117eB/117eA.hhr"
    arguments.add_argument('-hhr1',action='store',default='',help=hhrFile1Help)

    hhrFile2Help = "Optional. Full Path to the hhr file results for chain B. "\
            "Default checks for existence in inputDir/target/SPRING/ then "\
            "inputDir/target/ ie /home/user/117eA-117eB/SPRING/117eB.hhr or"\
            " /home/user/117eA-117eB/117eB.hhr"
    arguments.add_argument('-hhr2',action='store',default='',help=hhrFile2Help)

    hhsearchHelp = "Optional. If flag set spring will run hhsearch if results "\
            "file missing"
    arguments.add_argument('-hhsearch',action='store_true',
            default=springConfig.parameters['runHHsearch'],help=hhsearchHelp )

    arguments.add_argument('-W0',action="store",type=float,
            default=defaultParam["W0"],
            help='Weight for TMscore protion of Spring Score')
    arguments.add_argument('-W1',action="store",type=float,
            default=defaultParam["W1"],
            help='Weight for DFire portion of Spring Score')
    arguments.add_argument('-Norm',action="store",type=float,
            default=defaultParam["Normalize"],
            help='Value of normalization score')

    dockMonomerHelp = "Docks Top ranked monomer for each chain to to dimer "\
            "template framework."
    arguments.add_argument('-dockMono',action='store_true',default=False,help=dockMonomerHelp)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    return args 

def WriteList(queryList,fout_Path,batchSize):
    """
    For batch submission splits query List into a set of files for job
    submission.
    
    Returns:
        foutList ([str]): List of full paths of files.
    """
    qListOut = []
    count = 0
    for i in xrange(0,len(queryList), batchSize):
        foutName = fout_Path +'/springList_'+str(count+1)
        fout = open(foutName, 'w')
        outList = [query+'\n' for query in queryList[count*batchSize:(count+1)*batchSize] ]
        outStr = ''.join(outList)
        fout.write(outStr)
        fout.close()
        qListOut.append(foutName)
        count+=1
    return qListOut


if __name__ == "__main__":
    args = ParseCommandLine()
    query = args.q
    inputDirectory = args.iDir
    submitJob = args.job
    outputDirectory = args.oDir
    run_hhsearch = args.hhsearch
    W0 = args.W0
    W1 = args.W1
    Normalize = args.Norm
    #rstfile
    monfileA = args.c1
    monfileB = args.c2
    hhr1  = args.hhr1
    hhr2  = args.hhr2
    batchSize = springConfig.jobParameters['batchSize']
    dockMono = args.dockMono

    if batchSize < 1:
        batchSize = 1
    
    queryList = []
    queryStr = query
    if os.path.isfile(query):
        f = open(query)
        targetList = [line.strip() for line in f]
        f.close()
        if submitJob:
            if batchSize > 1:
                fout_Path = inputDirectory+'/record/'
                mkdir(fout_Path)
                queryList = WriteList(targetList, fout_Path, batchSize)
            else:
                queryList = targetList
        else:
            queryList.append(query)
    else:
        queryList.append(query)

    if len(queryList) > 1:
        monfileA = ''
        monfileB = ''
        hhr1  = ''
        hhr2  = ''
    
    for query in queryList:
        spring_object = SPRING(query,inputDirectory,
            outputDirectory=outputDirectory,
            W0=W0, W1=W1, Normalize=Normalize, fileModelChainA=monfileA,
            fileModelChainB=monfileB,hhrFileA=hhr1,hhrFileB=hhr2,
            runHHsearch = run_hhsearch, dockMonomer=dockMono)
        if submitJob:
            recordDir = ''
            if not os.path.isfile(queryStr):
                recordDir = inputDirectory+'/'+query+'/record/'
            else:
                recordDir = inputDirectory+'/record/'
            spring_object.Submit(recordDir)
        else:
            spring_object.RunSpring()
