from Bio import pairwise2
from os.path import isfile


class Alignment:
    def __init__(self, fileName):
        self.queryName = None
        self.queryStart = list()
        self.queryAlignment = list()
        self.templateName = None
        self.templateStart = list()
        self.templateAlignment = list()
        self.readFile(fileName)

    def readFile(self, fileName):
        with open(fileName) as file:
            for line in file:
                cols = line.split()
                if len(cols) > 1 and cols[0] == "Query":
                    self.queryName = cols[1].split()[0][0:14]
                if len(cols) > 1 and cols[0].startswith(">"):
                    self.templateName = cols[0][1:]
                if self.queryName and self.templateName:
                    if len(cols) > 2:
                        if cols[0] == "Q" and cols[1] == self.queryName:
                            self.queryStart.append(self.getStart(cols[2]))
                            self.queryAlignment.append(cols[3])
                        if cols[0] == "T" and cols[1] == self.templateName:
                            self.templateStart.append(self.getStart(cols[2]))
                            self.templateAlignment.append(cols[3])
                if len(cols) > 1 and cols[0] == "No" and cols[1] == "2":
                    break

    def createModel(self, templateChain):
        hhrMapping = self.mapSequence(templateChain)
        previousResidue = dict()
        for residueNumber in templateChain:
            templateResidue = templateChain[residueNumber]["residue"]
            previousResidue[residueNumber] = templateResidue
            templateChain[residueNumber]["residue"] = None
        tCount = 0
        for i in range(len(self.templateAlignment)):
            templateSequence = self.templateAlignment[i]
            querySequence = self.queryAlignment[i]
            queryStart = self.queryStart[i]
            n = len(querySequence)
            qCount = 0
            for j in range(n):
                qs = querySequence[j]
                ts = templateSequence[j]
                rs = hhrMapping[tCount]
                if rs != "-" and qs != "-" and ts != "-":
                    residueNumber = rs
                    if residueNumber in templateChain:
                        pr = previousResidue[residueNumber]
                        if pr != self.toThreeAmino(ts):
                            print("Warning: Ignoring mismatching residue [%s != %s]." % (pr, self.toThreeAmino(ts)))
                        templateChain[residueNumber]["residue"] = self.toThreeAmino(qs)
                        templateChain[residueNumber]["residueNumber"] = queryStart + qCount
                    else:
                        print("Warning: Skipping missing residue [%s]." % residueNumber)
                if qs != "-":
                    qCount = qCount + 1
                if ts != "-":
                    tCount = tCount + 1

    def mapSequence(self, templateChain):
        pdbSequence = ""
        for residueNumber in templateChain:
            templateResidue = templateChain[residueNumber]["residue"]
            pdbSequence = pdbSequence + self.toSingleAmino(templateResidue)
        hhrSequence = ""
        for i in range(len(self.templateAlignment)):
            templateSequence = self.templateAlignment[i]
            for s in templateSequence:
                if s != "-":
                    hhrSequence = hhrSequence + s
        alignments = pairwise2.align.globalxx(pdbSequence, hhrSequence)
        pdbAlignment = alignments[0].seqA
        hhrAlignment = alignments[0].seqB
        pCount = 0
        hhrMapping = []
        for i in range(len(pdbAlignment)):
            p = pdbAlignment[i]
            h = hhrAlignment[i]
            if h != "-":
                residueIndex = pCount if p != "-" else p 
                hhrMapping.append(residueIndex)
            if p != "-":
                pCount = pCount + 1
        sortedResidueIndex = sorted(templateChain.keys())
        for i in range(len(hhrMapping)):
            if hhrMapping[i] != "-":
                hhrMapping[i] = sortedResidueIndex[hhrMapping[i]]
        return hhrMapping

    def getStart(self, x):
        try:
            return int(x)
        except Exception:
            raise Exception("Invalid start index in alignment [%s]." % x)

    def toThreeAmino(self, seq):
        code = dict(G="GLY", A="ALA", V="VAL", L="LEU", I="ILE", M="MET", F="PHE", P="PRO", Y="TYR", W="TRP",
                    K="LYS", S="SER", C="CYS", N="ASN", Q="GLN", H="HIS", T="THR", E="GLU", D="ASP", R="ARG")
        return code[seq] if seq in code else "XXX"

    def toSingleAmino(self, seq):
        code = dict(GLY="G", ALA="A", VAL="V", LEU="L", ILE="I", MET="M", PHE="F", PRO="P", TYR="Y", TRP="W",
                    LYS="K", SER="S", CYS="C", ASN="N", GLN="Q", HIS="H", THR="T", GLU="E", ASP="D", ARG="R")
        return code[seq] if seq in code else "X"


def getTemplates(hhrFile, minScore):
    result = dict()
    topTemplate = None
    if isfile(hhrFile):
        with open(hhrFile) as file:
            for index, line in enumerate(file):
                if index > 8:
                    if not line.strip():
                        break
                    templateId = line[4:10]
                    templateScore = float(line[57:63])
                    if templateScore > minScore:
                        if topTemplate is None:
                            topTemplate = templateId
                        result[templateId] = templateScore
    return topTemplate, result


def getCrossReference(crossReferenceFile):
    crossReference = dict()
    with open(crossReferenceFile) as file:
        for line in file:
            columns = line.split()
            if len(columns) < 2:
                raise Exception("Invalid cross reference entry %s." % line)
            core = columns[0]
            partner = columns[-1]
            if core not in crossReference:
                crossReference[core] = []
            crossReference[core].append(partner)
    return crossReference
