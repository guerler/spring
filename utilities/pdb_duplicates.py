#! /usr/bin/env python3
import argparse


def formatEntry(line):
    return "%s.pdb" % line[:4].lower()


def main(args):
    existing = set()
    with open(args.a) as file:
        for line in file:
            entryA = formatEntry(line)
            existing.add(entryA)
    print("Found %s existing entries." % len(existing))
    entries = set()
    with open(args.b) as file:
        for line in file:
            entryB = formatEntry(line)
            if entryB not in existing:
                entries.add(entryB)
    print("Found %s new entries." % len(entries))
    with open("pdb_duplicates.log", "w") as file:
        for entry in entries:
            file.write("%s\n" % entry)
    file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Remove duplicate entries')
    parser.add_argument('-a', '--a', help='List of existing entries', required=True)
    parser.add_argument('-b', '--b', help='List of new entries', required=True)
    args = parser.parse_args()
    main(args)
