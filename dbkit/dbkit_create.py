#! /usr/bin/env python3
import argparse
from os import system
from os.path import getsize


def getIdentifiers(args):
    entries = set()
    with open(args.list) as file:
        for line in file:
            entry = line.split()[0]
            idLength = int(args.idlength)
            if idLength > 0:
                entry = entry[:idLength]
            if args.idcase == "lower":
                entry = entry.lower()
            elif args.idcase == "upper":
                entry = entry.upper()
            if args.idextension is not None:
                entry = "%s.%s" % (entry, args.idextension)
            entries.add(entry)
    return sorted(entries)


def main(args):
    entries = getIdentifiers(args)
    logFile = open(args.log, "w")
    logFile.write("Found %s entries.\n" % len(entries))
    outputIndex = args.index
    outputDatabase = args.database
    tempPath = args.temp.rstrip("/")
    tempFile = "%s/temp.pdb" % tempPath
    system("mkdir -p %s" % tempPath)
    system("rm -f %s" % outputDatabase)
    indexFile = open(outputIndex, 'w')
    start = 0
    for entryId in entries:
        logFile.write("Loading %s.\n" % entryId)
        system("wget -q -O %s %s%s" % (tempFile, args.url, entryId))
        tempSize = getsize(tempFile)
        if tempSize == 0:
            logFile.write("Entry `%s` not found.\n" % entryId)
        else:
            indexFile.write("%s\t%d\t%d\n" % (entryId, start, tempSize))
            start = start + tempSize + 1
            system("cat %s >> %s" % (tempFile, outputDatabase))
        logFile.flush()
    logFile.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='DBKit - Download and Merge files into a single file.')
    parser.add_argument('-l', '--list', help='List of entries', required=True)
    parser.add_argument('-u', '--url', help='Source Url', required=True)
    parser.add_argument('-t', '--temp', help='temp', required=True)
    parser.add_argument('-il', '--idlength', help='Format Identifier Length (integer)', required=False, default="0")
    parser.add_argument('-ic', '--idcase', help='Format Identifier Case (lower, upper)', required=False, default=None)
    parser.add_argument('-ie', '--idextension', help='Format Identifier Extension', required=False, default=None)
    parser.add_argument('-o', '--index', help='Output Database Index', required=True)
    parser.add_argument('-d', '--database', help='Output Database', required=True)
    parser.add_argument('-g', '--log', help="Log file", required=True)
    args = parser.parse_args()
    main(args)
