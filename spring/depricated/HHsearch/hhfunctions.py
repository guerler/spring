#!/usr/bin/env python
"""
Set of functions for handling hhsearch results.  Provides functions for making 
and writing to file hhsearch models.  If run as main function will write
hhsearch models to file.
"""
#python 2.x and 3.x cross compatibility
from __future__ import print_function, division
import sys
if sys.version_info[0] >= 3:
    xrange = range

import os
import pandas as pd
import argparse
import re
import math

import hhsearchConfig
rootDir = hhsearchConfig.parameters['rootDir']
sys.path.append(hhsearchConfig.parameters['pdbPath'])
from pdbchains import ReadFasta,PDBchains

class HHSEARCH_Functions:
    """
    This class is a module for extracting the data from hhesarch results file ie
    ('query.hhr') The module provides function for extracing the summary 
    information from the hhr files and making HHsearch Models.

    Attributes
        hhsearchResultsFile (str): Full path for hhsearch results file (.hhr)
        summary (pandas object):Stores summary information for each alignment in
            the hhsearchReultsFile.  More information can be found in the 
            comments of the ReadHHsearchSummary function
        alignmentData (dictionary): Contains information of each template to 
            query sequence alignment.  Such as the alignment string and 
            alignment positions.  More information can be found in comments
            of ReadHHsearch.
            alignmentsData (dictionary): Each value in the dictionary contains 4 
                            pieces of information
                QueryResidueNumberList (integer): list of residues in the 
                            alignment for the query
                QueryAlignedResidues   (string): string of aligned residues in 
                            the query sequence
                TemplateResidueNumberList (integer): list of residues in the 
                            alignment for the template
                TemplateAlignedResidues (string): string of aligned residues in 
                            the query sequence
            Example:
                alignmentData['12asA'] = [ [321,322,..,327], [E,S,V,P,S,L,L] , 
                                                [321,322,..,327], [E,S,V,P,S,L,L]]
    Methods
        ReadHHsearch:Extracts summary and alignmentData from hhsearchResultsFile
        ReadHHsearchSummary: Extracts summary data only form hhsearchResultsFile
        MakeModel: Make an HHsearch Model from the hhr file
        MakeModels: Make several HHsearch Models
        WriteModel: Write HHsearch Model to file.
    """
    
    def __init__(self,hhsearchResultsFile, GetAlignmentData = True):
        self.hhsearchResultsFile = hhsearchResultsFile
        self.summary = None
        self.alignmentData = None
        if GetAlignmentData:
            self.ReadHHsearch() #gets summary and alignment Data
        else:
            self.ReadHHsearchSummary() #gets summary only

    def WriteModel(self,templatePdbFile, templateName, outputName):
        """
        Writes HHsearch model to file.
        Arguments:
            templatePdbFile (str): Full path to template pdbfile
            templateName (str): Name of template.  ie if file is 12asA.pdb
                    templateName == 12asA
            outputName (str): Full path of output file
        """
        templateModel = self.MakeModel(templatePdbFile, templateName)
        templateModel.Write(outputName,write_template=True)

    def MakeModel(self,templatePdbFile,templateName):
        """
        Returns HHsearch Model. CA only
    
        Arguments:
            templatePdbFile (string): Full Path to template PDB file
            templateName (string): Name of template.  Ie if 12asA.pdb is pdb file
                12asA is the name of the template.  It must have same name
                in hhsearch database.
        Returns:
            HHsearchModel (pdbchains object): Template ca model.
        """
        if self.summary is None or self.alignmentData is None:
            self.ReadHHsearch()

        if templateName not in self.alignmentData:
            raise KeyError(templateName+' not in alignmentData check templateName or hhr file')
    
        template = PDBchains()
        template.Read(templatePdbFile, ca_only=True, allow_insert=False,
                        allow_alt_loc = False)
        template.Renumber()
        templateAlignmentData = self.alignmentData[templateName]
        QueryResidueNumberList, QueryAlignedSequence, TemplateResidueNumberList,TemplateAlignedSequence = templateAlignmentData
        #get subset of aligned residues in template
        templateModel = template.SliceResNum(0,TemplateResidueNumberList,
                                                ca_only = True)
        #switch residue number and sequence information.
        residueAbbreviationToName = templateModel.codeToabrev
        for i in xrange(0, len(templateModel.atom_info[0]) ):
            temp_res = templateModel.atom_info[0][i].res_name
            temp_num = templateModel.atom_info[0][i].res_num
            templateModel.atom_info[0][i].temp_res = temp_res
            templateModel.atom_info[0][i].temp_num = temp_num
            queryResName = residueAbbreviationToName[ QueryAlignedSequence[i] ]
            templateModel.atom_info[0][i].res_name = queryResName
            templateModel.atom_info[0][i].res_num = QueryResidueNumberList[i]
        return templateModel

    def MakeModels(self,templatePdbFiles, templateNames):
        """
        Returns a list of pdbchains objects containing hhsearch models
        Arguments:
            templatePdbFiles ([string]): List ofFull Path to template PDB file
            templateNames ([string]):List of Names of template.  Ie if 12asA.pdb
                is pdb file 12asA is the name of the template.  It must have 
                same name in hhsearch database.
        Note:
            The templateName and templatePdbFile should match. ie.
            if templatePdbFiles[0]=='12asA.pdb' then templateNames[0]=='12asA'
        Returns:
            HHsearchModels ([pdbchains object]): List of Template ca model in 
            pdbchains format
        """
        HHsearchModels = []
        for i in xrange(0,len(templateNames) ):
            templateFile = templatePdbFiles[i]
            templateName = templateNames[i]
            templateModel = self.MakeModel(templateFile, templateName)
            HHsearchModels.append(templateModel)
        return HHsearchModels

    def ReadHHsearch(self):
        """
        Reads hhsearch.hhr file and orgainzes alignment results.
    
        Arguments:
            hhsearchResultsFile (file): Full path of HHsearch results file
                for the query sequence.
            parseAlignment (bool): Return Alignment Information
        Returns (pandas DataFrame, Dictionary):
            summary (padas data frame): 
                Template (str): Name of template PDB
                Probability (float): Probability template to query alignment 
                    correct.
                Evalue (float): E-value for query to template alignment
                Score (float): Raw HHsearch score
                AlignendColumns (int): Number of aligned columns
                SeqID (float): Sequence identity of query to template
                    alignment. From 0 to 1.
                Similarity (float): HHsearch similarity score
                SumProbs (float): Another HHsearch score.
                SpringZscore (float): A log10 transformation of the E-value.
                    used in spring threading to represent quality
                    of threading result
        alignments (dictionary): sequence alignments from hhsearch.
        Example
            out = ReadHHsearch('query.hhr')  
        """
        f = open(self.hhsearchResultsFile)
        lines = f.readlines()
        query_name = lines[0].split()[1].strip() #query name
        #summary = {}
        alignments = {}
        rank = [] 
        current_template = ''
        is_hit = False  
        for i in xrange(0 ,len(lines)):
            if is_hit:
                if i+1 >= len(lines):
                    continue
                if lines[i+1][0] == '>':
                    is_hit = False
                    continue
                else:
                    continue

            if lines[i][0] == '>':
                current_template = lines[i][1:].strip()
                if current_template in alignments:
                    is_hit = True
                    continue
                rank.append(current_template)
                #summary[current_template] = lines[i+1]
                alignments[current_template] = [ [], [] ] #[ [query], [template] ]
                continue

            if current_template == '':
                continue
            
            #if lines[i][0] == 'Q' and query_name in lines[i]:  #if query_name in lines[i]: Append query alignment
            if lines[i][0] == 'Q':
                if 'ss_pred' not in lines[i] and 'Consensus' not in lines[i]:
                    alignments[current_template][0].append(lines[i])

            if current_template in lines[i] and lines[i][0] == 'T': #append template alingment
                alignments[current_template][1].append(lines[i])
        f.close()   
        self._ParseAlignments(alignments,rank) #modifies self.alignmentData
        self.ReadHHsearchSummary() #modifies self.summary

    def ReadHHsearchSummary(self):
        """
        Reads hhsearch.hhr file and obtains summary results.

        Modifies self.summary 
            (pandas DataFrame):
            Template (str): Name of template PDB
            Probability (float): Probability template to query alignment correct
            Evalue (float): E-value for query to template alignment
            Score (float): Raw HHsearch score
            AlignendColumns (int): Number of aligned columns
            SeqID (float): Sequence identity of query to template
                alignment. From 0 to 1.
            Similarity (float): HHsearch similarity score
            SumProbs (float): Another HHsearch score.
            SpringZscore (float): A log10 transformation of the E-value.
                used in spring threading to represent quality
                of threading result.
        Example
            hhfunctions = HHsearchFunctions('12asA.hhr')
            hhfunctions.ReadHHsearchSummary()    
        """
        f = open(self.hhsearchResultsFile)
        lines = f.readlines()
        summary = {}
        chainLength = {}
        rank = [] 
        current_template = ''
        for i in xrange(0,len(lines) ): 
            if lines[i][0] == '>':
                current_template = lines[i][1:].strip()
                if current_template in summary:
                    continue
                rank.append(current_template)
                summary[current_template] = lines[i+1]
                chainLength[current_template] = self._GetChainLength(lines,i+8) 
        #   return alignment, rank, summary
        f.close()
        summary = self._ParseSummary(summary,rank,chainLength)
        data = {'Template':rank,'Probability':summary[0],'Evalue':summary[1],
            'Score':summary[2],'AlignedColumns':summary[3],
            'SeqID':summary[4], 'Similarity':summary[5], 
            'SumProbs':summary[6],'SpringZscore':summary[7],
            'ChainLength':summary[8]}
        summary_out = pd.DataFrame(data,columns=['Template','Probability',
            'Evalue','Score','AlignedColumns',
            'SeqID','Similarity','SumProbs','SpringZscore','ChainLength'])
        self.summary = summary_out

    def _GetChainLength(self,lines,lineNum):
        """
        Extracts template chain length from alignment line.
        Arguments:
            lines ([str]): All lines in the hhr results file
            lineNum (int): Line number where template chain length is stored
            ie where last column is length:
            T 1il2A           139 KTRAKITSLVRRFMD-----DHGFLDIETPMLTKATPE------G-ARDYLVPSRVHKGKFYALPQSPQLFKQLLMMS-G  205 (585)
        Returns:
            chainLength (int): Length of template chain ie. 585 in the example
        """
        line = lines[lineNum]
        tmp = line.strip().split()
        chainLen = int( tmp[-1].replace('(','').replace(')','') )
        return chainLen

    def _ParseSummary(self,summary,rank,chainLength):
        """
        ReadHHsearch helper function.  Extracts the alignment 
        statistics from the hhsearch alignment
        """
        output = [[],[],[],[],[],[],[],[],[]]
        for template in rank:
            templateSummary = summary[template]
            chainLen = chainLength[template]
            #extract all numbers from templateSummary with the format below
            #'Probab=100.00  E-value=1.5e-62  Score=432.64  Aligned_cols=172  Identities=29%  Similarity=0.529  Sum_probs=164.7\n'
            values = re.findall(r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?',templateSummary)
            Probability = float(values[0])
            Evalue = float(values[1])
            Score = float(values[2])
            AlignedCols = float(values[3])
            PercentIdentities = float(values[4])/100
            Similarity = float(values[5])
            SumProbs = float(values[6])
            if Evalue == 0:
                Evalue = 7*(10**-107)
            SpringZscore = -1*math.log10(Evalue)
            output[0].append(Probability)
            output[1].append(Evalue)
            output[2].append(Score)
            output[3].append(AlignedCols)
            output[4].append(PercentIdentities)
            output[5].append(Similarity)
            output[6].append(SumProbs)
            output[7].append(SpringZscore)
            output[8].append(chainLen)
        return output

    def _ParseAlignments(self,alignments,rank):
        """
        Parses alignments dictionary and returns alignment structured alignment 
        information.

        Arguments:
            alignments (dictionary): Data structured described in Notes section.
            rank ([str]): Template list that is sorted based on similarity rank.
                ie rank[0] is the template that is most similar to the query
                sequence based on sequence similarity.

        Example Alignment segment from 12asA.hhr.
            Q ss_pred             HHhhhcC
            Q 12asA           321 ESVPSLL  327 (327)
            Q Consensus       321 ~~~~~~~  327 (327)
                          ++|||||
            T Consensus       321 ~~~~~~~  327 (327)
            T 11asA           321 ESVPSLL  327 (327)
            T ss_dssp             HHCCSCC
            T ss_pred             HHhhhcC
            Confidence            9999997

        Notes:
        The alignment dictionary will conatine each line containing sequence 
        information ie
        alignments['12asA'] == [ [ Q 11asA  321 ESVPSLL  327 (327) ], 
                                    [T 12asA  321 ESVPSLL  327 (327)] ]

        Returns: (Note: All list will be same length)
            alignmentsData (dictionary): Each value in the dictionary contains 4 
                            pieces of information
                QueryResidueNumberList (integer): list of residues in the 
                            alignment for the query
                QueryAlignedResidues   (string): string of aligned residues in 
                            the query sequence
                TemplateResidueNumberList (integer): list of residues in the 
                            alignment for the template
                TemplateAlignedResidues (string): string of aligned residues in 
                            the query sequence
            Example:
                alignmentData['12asA'] = [ [321,322,..,327], [E,S,V,P,S,L,L] , 
                                               [321,322,..,327],[E,S,V,P,S,L,L]]
        """
        alignmentData = {}
        rankNum = 0
        for template in alignments:
            queryData, templateData = alignments[template]
            QueryResidueNumberList = []
            QueryAlignedResidues = []
            TemplateResidueNumberList = []
            TemplateAlignedResidues = []
            for i in xrange(0,len(queryData) ):
                Qstart,QseqAlign,Qend = self._SplitAlignmentData(queryData[i])
                Tstart,TseqAlign,Tend =self._SplitAlignmentData(templateData[i])
                Qnum = Qstart
                Tnum = Tstart
                for j in xrange(0,len(QseqAlign) ):
                    Qpos = QseqAlign[j]
                    Tpos = TseqAlign[j]
                    if Qpos != '-' and Tpos != '-':
                        QueryResidueNumberList.append(Qnum)
                        QueryAlignedResidues.append(Qpos)
                        TemplateResidueNumberList.append(Tnum)
                        TemplateAlignedResidues.append(Tpos)
                        Qnum +=1
                        Tnum +=1
                    elif Tpos != '-':
                        Tnum += 1
                    else:
                        Qnum += 1
            alignmentData[template] = [QueryResidueNumberList, 
                            QueryAlignedResidues, TemplateResidueNumberList,
                            TemplateAlignedResidues]     
        self.alignmentData = alignmentData

    def _SplitAlignmentData(self, alignData):
        """
        Splits alignment data from hhsearch
        Arguments:
            alignData (string): ie, Q 12asA           321 ESVPSLL  327 (327)
        Returns:
            startAlign (int): Start of sequence alignment ie 321
            alignment (str): String of alignment ie ESVPSLL
            endAlign (int): End of sequence alignment ie 327
        """
        tmp = alignData.split()
        startAlign = int(tmp[2])
        alignment = tmp[3]
        endAlign = int(tmp[4])
        return startAlign, alignment, endAlign

def ParseCommandLine():
    """
    Parses command line for running hhfunctions.py
    Returns:
        args (Namespace for argparse): This namespace will contain all the 
                                       parameters from the command line

    """
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__),
                description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    arguments = parser.add_argument_group(title="Arguments")
    hhrFileHelp = "Full file path for hhr file generated by hhsearch"
    arguments.add_argument('-hhr',action='store',required=True,help=hhrFileHelp)
    templateFileHelp = "Full file path to template file"
    arguments.add_argument('-pdb', action='store',required=True,
            help=templateFileHelp)
    templateNameHelp = "Name of template pdb. If pdbfile is 12asA.pdb templateName is 12asA"
    arguments.add_argument('-tname',action='store',required=True,help=templateNameHelp)
    outputHelp = "Full path of output file"
    arguments.add_argument('-o',action='store',required=True,help=outputHelp)
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    return args

if __name__ == "__main__":

    args = ParseCommandLine()
    hhsearchResultsFile = args.hhr
    templateFile = args.pdb
    templateName = args.tname
    outFile = args.o
    hhfuncs = HHSEARCH_Functions(hhsearchResultsFile)
    hhfuncs.WriteModel(templateFile, templateName, outFile)
