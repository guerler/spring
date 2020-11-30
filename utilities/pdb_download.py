#! /usr/bin/env python3
import argparse
from os import system
from os.path import isfile, getsize


def main(args):
    pdbUrl = "https://files.rcsb.org/download/"
    pdbPath = args.pdbpath.rstrip("/")
    system("mkdir -p %s" % pdbPath)
    entries = set()
    with open(args.list) as file:
        for line in file:
            entries.add(line[:4].lower())
    print("Found %s template entries from `%s`." % (len(entries), args.list))
    for pdbId in sorted(entries):
        print("Loading %s..." % pdbId)
        pdbFile = "%s.pdb" % pdbId
        pdbPathFile = "%s/%s" % (pdbPath, pdbFile)
        if isfile(pdbPathFile):
            print("Skipping %s." % pdbId)
            continue
        system("wget -q -O %s %s%s" % (pdbPathFile, pdbUrl, pdbFile))
        if isfile(pdbPathFile) and getsize(pdbPathFile) == 0:
            print("Removing empty file %s." % pdbPathFile)
            system("rm %s" % pdbPathFile)
        print("Completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Downloads PDB files')
    parser.add_argument('-l', '--list', help='List of entries', required=True)
    parser.add_argument('-p', '--pdbpath', help='Path to files', required=True)
    args = parser.parse_args()
    main(args)
