#! /usr/bin/env python3
import argparse
from os import system
from os.path import getsize

from dbkit_package.DBKit import DBKit


def main(args):
    logFile = open(args.log, "w")
    outIndex = args.outindex
    outData = args.outdata
    if getsize(args.firstindex) > getsize(args.secondindex):
        firstIndex = args.firstindex
        firstData = args.firstdata
        secondIndex = args.secondindex
        secondData = args.seconddata
    else:
        firstIndex = args.secondindex
        firstData = args.seconddata
        secondIndex = args.firstindex
        secondData = args.firstdata
    system("cp %s %s" % (firstIndex, outIndex))
    system("cp %s %s" % (firstData, outData))
    firstEntries = set()
    with open(firstIndex, "r") as f:
        for line in f:
            name = line.split()[0]
            firstEntries.add(name)
    logFile.write("Detected %s entries.\n" % len(firstEntries))
    secondEntries = list()
    with open(secondIndex, "r") as f:
        for line in f:
            name = line.split()[0]
            secondEntries.append(name)
    tempFile = "temp.dat"
    count = 0
    dbkit = DBKit(secondIndex, secondData)
    for secondKey in secondEntries:
        if secondKey not in firstEntries:
            dbkit.createFile(secondKey, tempFile)
            entrySize = getsize(tempFile)
            currentSize = getsize(outData)
            system("cat %s >> %s" % (tempFile, outData))
            system("echo '%s\t%s\t%s' >> %s" % (secondKey, currentSize, entrySize, outIndex))
            count = count + 1
        else:
            logFile.write("Skipping existing entry %s.\n" % secondKey)
    logFile.write("Added %s entries.\n" % count)
    logFile.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='DBKit - Merge database pair.')
    parser.add_argument('-fi', '--firstindex', help='First Index file', required=True)
    parser.add_argument('-fd', '--firstdata', help='First Data file', required=True)
    parser.add_argument('-si', '--secondindex', help='Second Index file', required=True)
    parser.add_argument('-sd', '--seconddata', help='Second Data file', required=True)
    parser.add_argument('-oi', '--outindex', help='Output Index file', required=True)
    parser.add_argument('-od', '--outdata', help='Output Data file', required=True)
    parser.add_argument('-log', '--log', help='Log file', required=True)
    args = parser.parse_args()
    main(args)
