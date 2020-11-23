#! /usr/bin/env python3
import argparse
import os


def main(args):
    pdbUrl = "https://files.rcsb.org/download/"
    pdbPath = args.pdbpath.rstrip("/")
    os.system("mkdir -p %s" % pdbPath)
    entries = list()
    with open(args.list) as file:
        for line in file:
            entries.append(line.strip())
    print("Found %s template entries from `%s`." % (len(entries), args.list))

    for entryId in entries:
        pdbId = entryId[:4].lower()
        print("Loading %s..." % pdbId)
        pdbFile = "%s.pdb" % pdbId
        pdbPathFile = "%s/%s" % (pdbPath, pdbFile)
        if os.path.isfile(pdbPathFile):
            print("Skipping %s." % pdbId)
            continue
        os.system("wget -O %s %s%s" % (pdbFile, pdbUrl, pdbFile))
        os.system("mv %s %s" % (pdbFile, pdbPathFile))
        print("Completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Downloads PDB files')
    parser.add_argument('-l', '--list', help='List of entries', required=True)
    parser.add_argument('-p', '--pdbpath', help='Path to files', required=True)
    args = parser.parse_args()
    main(args)