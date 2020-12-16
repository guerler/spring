#! /usr/bin/env python3
import argparse
from os import system
from os.path import getsize

from dbkit_package.DBKit import DBKit


def main(args):
    logFile = open(args.log, "w")
    outIndex = args.outindex
    outData = args.outdata
    entries = list()
    with open(args.list, "r") as f:
        for line in f:
            name = line.split()[0]
            entries.append(name)
    logFile.write("Detected %s entries.\n" % len(entries))
    tempFile = "temp.dat"
    count = 0
    dbkit = DBKit(args.index, args.database)
    for entry in sorted(entries):
        success = dbkit.createFile(entry, tempFile)
        if success:
            currentSize = getsize(outData)
            entrySize = getsize(tempFile)
            system("cat %s >> %s" % (tempFile, outData))
            system("echo '%s\t%s\t%s' >> %s" % (entry, currentSize, entrySize, outIndex))
            count = count + 1
        else:
            logFile.write("Entry %s not found.\n" % entry)
    logFile.write("Extracted %s entries.\n" % count)
    logFile.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='DBKit - Merge database pair.')
    parser.add_argument('-l', '--list', help='List of entries to be extracted', required=True)
    parser.add_argument('-i', '--index', help='Database Index file (ffindex)', required=True)
    parser.add_argument('-d', '--database', help='Database Data file (ffdata)', required=True)
    parser.add_argument('-oi', '--outindex', help='Output Index file', required=True)
    parser.add_argument('-od', '--outdata', help='Output Data file', required=True)
    parser.add_argument('-g', '--log', help='Log file', required=True)
    args = parser.parse_args()
    main(args)
