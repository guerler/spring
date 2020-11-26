#! /usr/bin/env python3
import argparse
import os


def main(args):
    hhrPath = args.hhr_path.rstrip("/")
    outPath = args.output_path.rstrip("/")
    os.system("mkdir -p %s" % outPath)
    with open(args.input, "r") as file:
        for line in file:
            param = line.split()
            a = "%s/%s" % (hhrPath, param[0])
            b = "%s/%s" % (hhrPath, param[1])
            o = "%s/%s.%s" % (outPath, param[0], param[1])
            cmdString = "./spring_model.py -a '%s' -b '%s' -c '%s' -i '%s' -d '%s' -o '%s'" % (a, b, args.cross, args.index, args.database, o)
            os.system(cmdString)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a 3D model from HH-search results.')
    parser.add_argument('-p', '--pairs', help='Result table from min-Z evaluation', required=True)
    parser.add_argument('-hp', '--hhr_path', help='Path to HHR files', required=True)
    parser.add_argument('-op', '--output_path', help='Output directory', required=True)
    parser.add_argument('-i', '--index', help='PDB Database Index file (ffindex)', required=True)
    parser.add_argument('-d', '--database', help='PDB Database files (ffdata)', required=True)
    parser.add_argument('-c', '--cross', help='PDB Cross Reference', required=True)
    args = parser.parse_args()
    main(args)
