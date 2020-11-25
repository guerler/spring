#! /usr/bin/env python3
import argparse
import os


from spring_package.Molecule import Molecule


def getId(line):
    line = line.strip()
    return line[:4].upper() + line[4:6]


def getPDB(pdbPath, line):
    pdb = line[:4].lower()
    pdbChain = line[5:6]
    pdbFile = "%s/%s.pdb" % (pdbPath, pdb)
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


def getHomologue(queryFile, queryResultFile, databaseFile):
    if not os.path.isfile(queryResultFile):
        os.system("psiblast -query %s -db %s -out %s" % (queryFile,
                  databaseFile, queryResultFile))
    maxMatch = None
    try:
        with open(queryResultFile) as file:
            for i in range(38):
                line = next(file)
            maxMatch = getId(line.split()[0])
    except Exception:
        return None
    return maxMatch


def main(args):
    pdbPath = args.pdbpath.rstrip("/")
    temp = args.temp.rstrip("/")
    templates = set()
    os.system("mkdir -p %s" % temp)
    templateSequenceFile = "%s/templates.fasta" % temp
    if not os.path.isfile(templateSequenceFile):
        templateSequences = open(templateSequenceFile, "w")
        with open(args.list) as file:
            for rawId in file:
                templateId = getId(rawId)
                templates.add(templateId)
                pdbFile, pdbChain = getPDB(pdbPath, templateId)
                try:
                    templateMol = Molecule(pdbFile)
                    templateSeq = templateMol.getSequence(pdbChain)
                    templateSequences.write(">%s\n" % templateId)
                    templateSequences.write("%s\n" % templateSeq)
                except Exception:
                    print("Warning: File not found [%s]." % pdbFile)
        templateSequences.close()
        os.system("makeblastdb -in %s -dbtype prot" % templateSequenceFile)
    else:
        print("Using existing sequences for templates [%s]." %
              templateSequenceFile)
    print("Found %s template entries from `%s`." % (len(templates), args.list))

    crossReference = list()
    with open(args.cross) as file:
        for line in file:
            cols = line.split()
            if len(cols) != 2:
                raise Exception("Invalid line in crossreference [%s]." % line)
            crossReference.append(dict(core=cols[0], partner=cols[1]))
    print("Loaded crossreference with %d entries." % len(crossReference))

    for refEntry in crossReference:
        coreId = refEntry["core"]
        partnerId = refEntry["partner"]
        partnerFile = "%s/%s.fasta" % (temp, partnerId)
        partnerResultFile = "%s.result" % partnerFile
        if partnerId in templates:
            refEntry["match"] = partnerId
            print("Found partner in template list alignment [%s]" % partnerId)
        else:
            print("Processing %s." % partnerId)
            if not os.path.isfile(partnerResultFile):
                pdbFile, pdbChain = getPDB(pdbPath, partnerId)
                partnerMol = Molecule(pdbFile)
                partnerSeq = partnerMol.getSequence(pdbChain)
                with open(partnerFile, "w") as partnerFasta:
                    partnerFasta.write(">%s\n" % partnerId)
                    partnerFasta.write("%s" % partnerSeq)
            else:
                print("Using existing results. [%s]" % partnerId)
            matchedId = getHomologue(partnerFile, partnerResultFile,
                                     templateSequenceFile)
            if matchedId is None:
                print("Warning: Failed alignment [%s]" % partnerId)
            else:
                print("Found matching entry %s." % matchedId)
                refEntry["match"] = matchedId

    finalSet = set()
    for refEntry in crossReference:
        coreId = refEntry["core"]
        partnerId = refEntry["partner"]
        if "match" in refEntry:
            entry = "%s\t%s" % (coreId, refEntry["match"])
            if entry not in finalSet:
                finalSet.add(entry)
        else:
            print("Warning: Skipping failed missing partner match [%s, %s]." % (coreId, partnerId))
    print("Found %s cross reference entries." % len(finalSet))

    with open(args.output, 'w') as output_file:
        for entry in sorted(finalSet):
            output_file.write("%s\n" % entry)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Maps binding partners to template library')
    parser.add_argument('-l', '--list', help='List of template entries [PDB_CHAIN]', required=True)
    parser.add_argument('-c', '--cross', help='Cross reference (unmapped)', required=True)
    parser.add_argument('-p', '--pdbpath', help='Path to PDB files [PDB.pdb]', required=True)
    parser.add_argument('-o', '--output', help='Cross reference', required=True)
    parser.add_argument('-t', '--temp', help='Temporary directory', required=False, default="temp/")
    args = parser.parse_args()
    main(args)
