#! /usr/bin/env python3
import argparse
import os

def getIdentifier(line):
    pdbChain = line.split("_")
    if len(pdbChain) != 2:
        raise Exception("Invalid monomer identifier [%s]" % pdbChain)
    pdb = pdbChain[0]
    chain = pdbChain[1]
    return pdb, chain
 
def main(args):
    pdbPath = args.pdb_path.rstrip("/")
    hhrPath = args.hhr_path.rstrip("/")
    outPath = args.output_path.rstrip("/")
    os.system("mkdir -p %s" % outPath)
    with open(args.input, "r") as file:
        for line in file:
            param = line.split()
            ar = "%s/%s%s" % (hhrPath, param[0], ".txt")
            br = "%s/%s%s" % (hhrPath, param[1], ".txt")
            at, ac = getIdentifier(param[3])
            bt, bc = getIdentifier(param[4])
            ct, cc = getIdentifier(param[5])
            at = "%s/%s%s" % (pdbPath, at, ".pdb")
            bt = "%s/%s%s" % (pdbPath, bt, ".pdb")
            ct = "%s/%s%s" % (pdbPath, ct, ".pdb")
            o = "%s/%s.%s" % (outPath, param[0], param[1])
            cmdString = "./spring_model.py -ar '%s' -at '%s' -ac '%s' -br '%s' -bt '%s' -bc '%s' -ct '%s' -cc '%s' -o '%s'" % (ar, at, ac, br, bt, bc, ct, cc, o)
            os.system(cmdString)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a 3D model from HH-search results.')
    parser.add_argument('-i', '--input', help='Result table from min-Z evaluation', required=True)
    parser.add_argument('-hp', '--hhr_path', help='Path to HHR files', required=True)
    parser.add_argument('-pp', '--pdb_path', help='Path to PDB files', required=True)
    parser.add_argument('-tp', '--temp', help='Temporary directory', required=False, default="temp/")
    parser.add_argument('-op', '--output_path', help='Output directory', required=True)
    args = parser.parse_args()
    main(args)
 