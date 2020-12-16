#! /usr/bin/env python3
import argparse
from os import mkdir
from os.path import isdir

from spring_package.DBKit import DBKit
from spring_package.Modeller import createModel


class ModelArguments:
    def __init__(self, args):
        self.log = args.log
        self.index = args.index
        self.database = args.database
        self.cross = args.cross

    def set(self, a, b, output):
        self.a = a
        self.b = b
        self.output = output


def main(args):
    modelArgs = ModelArguments()
    outPath = args.output_path.rstrip("/")
    if not isdir(outPath):
        mkdir(outPath)
    if not isdir("temp"):
        mkdir("temp")
    aFile = "temp/a.hhr"
    bFile = "temp/b.hhr"
    dbkit = DBKit(args.hhr_index, args.hhr_database)
    with open(args.pairs, "r") as file:
        for line in file:
            param = line.split()
            aIdentifier = param[0]
            bIdentifier = param[1]
            success = dbkit.createFile(aIdentifier, aFile)
            if not success:
                print("Failed to retrieve entry %s." % aIdentifier)
                continue
            success = dbkit.createFile(bIdentifier, bFile)
            if not success:
                print("Failed to retrieve entry %s." % bIdentifier)
                continue
            output = "%s/%s.%s.pdb" % (outPath, aIdentifier, bIdentifier)
            modelArgs.set(a=aFile, b=bFile, output=output)
            createModel(modelArgs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create 3D models from HH-search results.')
    parser.add_argument('-p', '--pairs', help='Interaction table e.g. from min-Z evaluation (2-columns)', required=True)
    parser.add_argument('-ih', '--hhr_index', help='HHR Index database file (ffindex)', required=True)
    parser.add_argument('-dh', '--hhr_database', help='HHR Database file (ffdata)', required=True)
    parser.add_argument('-i', '--index', help='PDB Database Index file (ffindex)', required=True)
    parser.add_argument('-d', '--database', help='PDB Database file (ffdata)', required=True)
    parser.add_argument('-c', '--cross', help='PDB Cross Reference', required=True)
    parser.add_argument('-g', '--log', help='Log file', required=True)
    args = parser.parse_args()
    main(args)
