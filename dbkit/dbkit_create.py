#! /usr/bin/env python3
import argparse
from os import listdir, system
from os.path import isfile, join, getsize


def createIndex(inputPath, outputIndex, outputDatabase):
    names = sorted([f for f in listdir(inputPath) if isfile(join(inputPath, f))])
    files = [join(args.path, name) for name in names]
    sizes = [getsize(f) for f in files]
    start = 0
    system("rm -f %s" % outputDatabase)
    with open(outputIndex, 'w') as output_file:
        for i in range(len(names)):
            name = names[i]
            size = sizes[i]
            file = files[i]
            end = start + size
            output_file.write("%s\t%d\t%d\n" % (name, start, end))
            start = end + 1
            system("cat %s >> %s" % (file, outputDatabase))


def downloadFiles(entries, sourceUrl, outputPath, logFile):
    outputPath = outputPath.rstrip("/")
    system("mkdir -p %s" % outputPath)
    for entryId in entries:
        logFile.write("Loading %s.\n" % entryId)
        fileName = "%s/%s" % (outputPath, entryId)
        if isfile(fileName):
            logFile.write("Skipping %s.\n" % entryId)
            continue
        system("wget -q -O %s %s%s" % (fileName, sourceUrl, entryId))
        if getsize(fileName) == 0:
            logFile.write("Entry not found.\n")
            system("rm %s" % fileName)
        else:
            logFile.write("Completed.\n")


def main(args):
    logFile = open(args.log, "w")
    entries = list()
    with open(args.list) as file:
        for line in file:
            entry = line.split()[0]
            if args.idlength is not None:
                entry = entry[:args.idlength]
            if args.idcase == "lower":
                entry = entry.lower()
            elif args.idcase == "upper":
                entry = entry.upper()
            if args.idextension is not None:
                entry = "%s.%s" % (entry, args.idextension)
            entries.append(entry)
    logFile.write("Found %s entries from `%s`.\n" % (len(entries), entries))
    downloadFiles(entries, args.url, args.path, logFile=logFile)
    createIndex(args.path, args.index, args.database)
    logFile.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='DBKit - Download and Merge files into a single file.')
    parser.add_argument('-l', '--list', help='List of entries', required=True)
    parser.add_argument('-u', '--url', help='Source Url', required=True)
    parser.add_argument('-p', '--path', help='Path to files', required=True)
    parser.add_argument('-x', '--index', help='Output Database Index', required=True)
    parser.add_argument('-d', '--database', help='Output Database', required=True)
    parser.add_argument('-il', '--idlength', help='Format Identifier Length (integer)', type=int, required=False, default=None)
    parser.add_argument('-ic', '--idcase', help='Format Identifier Case (lower, upper)', required=False, default=None)
    parser.add_argument('-ie', '--idextension', help='Format Identifier Extension', required=False, default=None)
    parser.add_argument('-o', '--log', help="Log file", required=True)
    args = parser.parse_args()
    main(args)
