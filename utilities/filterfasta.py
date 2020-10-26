#! /usr/bin/env python3
import argparse
import os

def main(args):
	names = set()
	with open(args.list) as file:
		for index, line in enumerate(file):
			names.add(line.strip())
	print ("Loaded %s names from `%s`." % (len(names), args.list))
	with open(args.output, 'w') as output_file:
		with open(args.fasta) as file:
			for index, line in enumerate(file):
				if line.startswith(">"):
					name = line.split()[0][1:]
					if name in names:
						output_file.write("%s" % line)
						nextLine = next(file)
						while nextLine and not nextLine.startswith(">"):
							output_file.write("%s" % nextLine)
							nextLine = next(file)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='This script filters sequences by identifier from a fasta file.')
	parser.add_argument('-l', '--list', help='Fasta containing mulitple sequences', required=True)
	parser.add_argument('-f', '--fasta', help='Fasta input file', required=True)
	parser.add_argument('-o', '--output', help='Output file containing filtered sequences', required=True)
	args = parser.parse_args()
	main(args)