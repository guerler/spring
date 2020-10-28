#! /usr/bin/env python3
import argparse
import os

def main(args):
	names = set()
	with open(args.list) as file:
		for line in file:
			names.add(line[0:args.idlength].upper())
	print ("Loaded %s names from `%s`." % (len(names), args.list))
	with open(args.output, 'w') as output_file:
		with open(args.fasta) as file:
			nextLine = next(file, None)
			while nextLine:
				if nextLine.startswith(">"):
					name = nextLine.split()[0][1:args.idlength+1].upper()
					if name in names:
						output_file.write(">%s\n" % name)
						nextLine = next(file)
						while nextLine and not nextLine.startswith(">"):
							output_file.write("%s" % nextLine)
							nextLine = next(file, None)
						names.remove(name)
					else:
						nextLine = next(file, None)
				else:
					nextLine = next(file, None)
	if len(names) > 0:
		print("Warning: Missing %s sequences." % len(names))
		print(", ".join(names))

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='This script filters sequences by identifier from a fasta file.')
	parser.add_argument('-l', '--list', help='List of entries', required=True)
	parser.add_argument('-f', '--fasta', help='Fasta input file', required=True)
	parser.add_argument('-o', '--output', help='Output file containing filtered sequences', required=True)
	parser.add_argument('-idx', '--idlength', help='Length of identifer', type=int, default=6)
	args = parser.parse_args()
	main(args)