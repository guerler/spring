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
    os.system("mkdir -p %s" % args.temp)
    with open(args.input, "r") as file:
        pdbPath = args.pdb_path.rstrip("/")
        hhrPath = args.hhr_path.rstrip("/")
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
            cmdString = "./spring_model.py -ar '%s' -at '%s' -ac '%s' -br '%s' -bt '%s' -bc '%s' -ct '%s' -cc '%s'" % (ar, at, ac, br, bt, bc, ct, cc)
            print(cmdString)
            os.system(cmdString)
            return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a 3D model from HH-search results.')
    parser.add_argument('-i', '--input', help='Result table from min-Z evaluation', required=True)
    parser.add_argument('-hp', '--hhr_path', help='Path to HHR files', required=True)
    parser.add_argument('-pp', '--pdb_path', help='Path to PDB files', required=True)
    parser.add_argument('-tp', '--temp', help='Temporary directory', required=False, default="temp/")
    parser.add_argument('-o', '--output', help='Output directory', required=False)
    args = parser.parse_args()
    main(args)
 