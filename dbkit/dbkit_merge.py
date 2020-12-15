#! /usr/bin/env python3
import argparse
from os import system
from os.path import getsize


def createFile(sourceData, start, size, outputName):
    with open(sourceData) as file:
        file.seek(int(start))
        content = file.read(int(size))
        outputFile = open(outputName, "w")
        outputFile.write(content)
        outputFile.close()


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
    secondEntries = dict()
    with open(secondIndex, "r") as f:
        for line in f:
            lineSplit = line.split()
            if len(lineSplit) < 3:
                raise Exception("Invalid Index entry found: %s." % line)
            name = lineSplit[0]
            start = lineSplit[1]
            size = lineSplit[2]
            secondEntries[name] = [start, size]
    tempFile = "temp.dat"
    count = 0
    for secondKey in secondEntries:
        if secondKey not in firstEntries:
            secondLocation = secondEntries[secondKey]
            createFile(secondData, secondLocation[0], secondLocation[1], tempFile)
            currentSize = getsize(outData)
            system("cat %s >> %s" % (tempFile, outData))
            system("echo '%s\t%s\t%s' >> %s" % (secondKey, currentSize, secondLocation[1], outIndex))
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
