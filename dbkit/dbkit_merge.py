#! /usr/bin/env python3
import argparse
from os import system
from os.path import getsize


def createFile(self, sourceData, start, size, outputName):
    with open(sourceData) as file:
        file.seek(start)
        content = file.read(size)
        outputFile = open(outputName, "w")
        outputFile.write(content)
        outputFile.close()


def main(args):
    logFile = open(args.log, "w")
    outIndex = args.outindex
    outData = args.outData
    if getsize(args.firstindex) > getsize(args.secondindex):
        firstIndex = args.firstindex
        firstData = args.firstData
        secondIndex = args.secondindex
        secondData = args.seconddata
    else:
        firstIndex = args.secondindex
        firstData = args.secondData
        secondIndex = args.firstindex
        secondData = args.firstdata
    system("cp %s %s" % (firstIndex, outIndex))
    system("cp %s %s" % (firstData, outData))
    firstEntries = set()
    with open(firstIndex, "r") as f:
        for line in f:
            name = line.split[0]
            firstEntries.add(name)
    secondEntries = dict()
    with open(secondIndex, "r") as f:
        for line in f:
            lineSplit = line.split()
            name = lineSplit[0]
            start = lineSplit[1]
            size = lineSplit[2]
            secondEntries[name] = [start, size]
    tempFile = "temp.dat"
    count = 0
    for secondKey in secondEntries:
        if secondKey not in firstEntries:
            secondLocation = secondEntries[secondKey]
            createFile(secondData, secondLocation["start"], secondLocation["size"], tempFile)
            currentSize = getsize(outData)
            system("cat %s > %s" % (tempFile, outData))
            system("cat '%s\t%s\n' > %s" % (currentSize, secondLocation["size"], outIndex))
            count = count + 1
    logFile.write("Merged %s entries." % count)    
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
