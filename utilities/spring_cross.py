#! /usr/bin/env python3
import argparse
from spring_package.Molecule import Molecule


def getId(line):
    line = line.split()[0]
    return line[:4].upper() + line[4:6]


def createFile(identifier, databaseIndex, database, outputName):
    start = -1
    end = -1
    with open(databaseIndex) as file:
        for line in file:
            cols = line.split()
            if identifier == cols[0]:
                start = cols[1]
                end = cols[2]
                break
    if start != -1:
        with open(database) as file:
            file.seek(start)
            content = file.read(end - start)
            outputFile = open(outputName, "w")
            outputFile.write(content)
            outputFile.close()
        return True
    else:
        return False


def main(args):
    pdbCount = 0
    pdbDatabase = args.pdbpath
    partnerList = set()
    entries = list()
    with open(args.list) as file:
        for line in file:
            entries.append(getId(line))
    print("Found %s template entries from `%s`." % (len(entries), args.list))
    for entryId in entries:
        pdb = "%s.pdb" % entryId[:4].lower()
        pdbChain = entryId[5:6]
        pdbFile = "%s/%s" % (args.temp, pdb)
        createFile(pdb, args.pdbindex, args.pdbpath, pdbFile)
        try:
            mol = Molecule(pdbFile)
        except Exception:
            print("Warning: File '%s' not found" % pdbFile)
            continue
        pdbCount = pdbCount + 1
        print("Processing %s, chain %s." % (pdbFile, pdbChain))
        print("Found %d biomolecule(s)." % len(mol.biomol.keys()))
        for biomolNumber in mol.biomol:
            if biomolNumber == 0:
                print("Processing biomolecule.")
                bioMolecule = mol
            else:
                print("Processing biomolecule %d." % biomolNumber)
                bioMolecule = mol.createUnit(biomolNumber)
            nChains = len(bioMolecule.calpha.keys())
            print("Found %d chain(s)." % nChains)
            if nChains > 1 and pdbChain in bioMolecule.calpha:
                for bioChain in bioMolecule.calpha:
                    if bioChain == pdbChain:
                        continue
                    partnerPdbChain = "%s_%s" % (pdb.upper(), bioChain[:1])
                    partnerList.add("%s\t%s" % (entryId, partnerPdbChain))
            else:
                print("Skipping: Chain not found or single chain [%s]." %
                      pdbChain)
    with open(args.output, 'w') as output_file:
        for entry in sorted(partnerList):
            output_file.write("%s\n" % entry)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='List filtering.')
    parser.add_argument('-l', '--list', help='List of PDB chains [PDB_CHAIN]', required=True)
    parser.add_argument('-i', '--index', help='PDB Database Index file (dbkit_index)', required=True)
    parser.add_argument('-d', '--database', help='PDB Database files (dbkit)', required=True)
    parser.add_argument('-o', '--output', help='Output file', required=True)
    parser.add_argument('-t', '--temp', help='Temporary Directory', required=True))
    args = parser.parse_args()
    main(args)
