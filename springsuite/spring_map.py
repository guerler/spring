#! /usr/bin/env python3
import argparse
from os import system
from os.path import isfile

from spring_package.DBKit import DBKit
from spring_package.Molecule import Molecule
from spring_package.Utilities import getId, getChain, getName


def getPDB(line, pdbDatabase):
    pdb = getName(line)
    pdbChain = getChain(line)
    pdbFile = "%s/temp.pdb" % args.temp
    pdbDatabaseId = "%s.pdb" % pdb
    pdbDatabase.createFile(pdbDatabaseId, pdbFile)
    return pdbFile, pdbChain


def getSequences(fileName):
    sequences = dict()
    with open(fileName) as file:
        for line in file:
            if line.startswith(">"):
                name = getId(line.split()[0][1:])
                nextLine = next(file)
                sequences[name] = nextLine
    return sequences


def findMatch(identifier, databaseFile, pdbDatabase):
    fastaFile = "temp/%s.fasta" % identifier
    resultFile = "%s.result" % fastaFile
    if not isfile(resultFile):
        pdbFile, pdbChain = getPDB(identifier, pdbDatabase)
        mol = Molecule(pdbFile)
        seq = mol.getSequence(pdbChain)
        with open(fastaFile, "w") as fasta:
            fasta.write(">%s\n" % identifier)
            fasta.write("%s" % seq)
        system("psiblast -query %s -db %s -out %s" % (fastaFile, databaseFile, resultFile))
    maxMatch = None
    try:
        with open(resultFile) as file:
            for i in range(38):
                line = next(file)
            maxMatch = getId(line.split()[0])
    except Exception:
        return None
    return maxMatch


def main(args):
    logFile = open(args.log, "w")
    temp = args.temp.rstrip("/")
    templates = set()
    system("mkdir -p %s" % temp)
    templateSequenceFile = "%s/templates.fasta" % temp
    pdbDatabase = DBKit(args.index, args.database)
    if not isfile(templateSequenceFile):
        templateSequences = open(templateSequenceFile, "w")
        with open(args.list) as file:
            for rawId in file:
                templateId = getId(rawId)
                templates.add(templateId)
                pdbFile, pdbChain = getPDB(templateId, pdbDatabase)
                try:
                    templateMol = Molecule(pdbFile)
                    templateSeq = templateMol.getSequence(pdbChain)
                    templateSequences.write(">%s\n" % templateId)
                    templateSequences.write("%s\n" % templateSeq)
                except Exception:
                    logFile.write("Warning: File not found [%s].\n" % pdbFile)
        templateSequences.close()
        system("makeblastdb -in %s -dbtype prot" % templateSequenceFile)
    else:
        logFile.write("Using existing sequences for templates [%s].\n" % templateSequenceFile)
    logFile.write("Found %s template entries from `%s`.\n" % (len(templates), args.list))
    logFile.flush()

    crossReference = list()
    with open(args.cross) as file:
        for line in file:
            cols = line.split()
            if len(cols) != 2:
                raise Exception("Invalid line in crossreference [%s]." % line)
            crossReference.append(dict(core=cols[0], partner=cols[1]))
    logFile.write("Loaded crossreference with %d entries.\n" % len(crossReference))
    logFile.flush()

    for refEntry in crossReference:
        coreId = refEntry["core"]
        partnerId = refEntry["partner"]
        logFile.write("Processing %s.\n" % partnerId)
        if partnerId in templates:
            partnerMatch = partnerId
        else:
            partnerMatch = findMatch(partnerId, templateSequenceFile, pdbDatabase)
        if partnerMatch is None:
            logFile.write("Warning: Failed alignment [%s].\n" % partnerId)
        else:
            logFile.write("Found matching entry %s.\n" % partnerMatch)
            refEntry["partnerMatch"] = partnerMatch
        logFile.flush()

    finalSet = set()
    for refEntry in crossReference:
        coreId = refEntry["core"]
        partnerId = refEntry["partner"]
        if "partnerMatch" in refEntry:
            entry = "%s\t%s" % (coreId, refEntry["partnerMatch"])
            if entry not in finalSet:
                finalSet.add(entry)
        else:
            logFile.write("Warning: Skipping failed missing partner match [%s, %s].\n" % (coreId, partnerId))
    logFile.write("Found %s cross reference entries.\n" % len(finalSet))
    logFile.close()

    with open(args.output, 'w') as output_file:
        for entry in sorted(finalSet):
            output_file.write("%s\n" % entry)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Maps binding partners to template library')
    parser.add_argument('-l', '--list', help='List of template entries `PDB_CHAIN`', required=True)
    parser.add_argument('-i', '--index', help='PDB Database Index file (ffindex)', required=True)
    parser.add_argument('-d', '--database', help='PDB Database files (ffdata)', required=True)
    parser.add_argument('-c', '--cross', help='Cross reference (unmapped)', required=True)
    parser.add_argument('-o', '--output', help='Cross reference', required=True)
    parser.add_argument('-g', '--log', help='Log File', required=True)
    parser.add_argument('-t', '--temp', help='Temporary directory', required=True)
    args = parser.parse_args()
    main(args)
