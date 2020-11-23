#! /usr/bin/env python3
import argparse
from os import listdir, system
from os.path import isfile, join, getsize


def createIndex(inputPath, outputIndex, outputDatabase):
    names = sorted([f for f in listdir(inputPath) if isfile(join(inputPath, f))])
    files = [join(args.path, name) for name in names]
    sizes = [getsize(f) for f in files]
    start = 0
    system("rm %s" % outputDatabase)
    with open(outputIndex, 'w') as output_file:
        for i in range(len(names)):
            name = names[i]
            size = sizes[i]
            file = files[i]
            end = start + size
            output_file.write("%s\t%d\t%d\n" % (name, start, end))
            start = end + 1
            system("cat %s >> %s" % (file, outputDatabase))


def main(args):
    createIndex(args.path, args.index, args.database)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merges files.')
    parser.add_argument('-p', '--path', help='Path to files', required=True)
    parser.add_argument('-x', '--index', help='Output Database Index', required=True)
    parser.add_argument('-d', '--database', help='Output Database', required=True)
    args = parser.parse_args()
    main(args)
