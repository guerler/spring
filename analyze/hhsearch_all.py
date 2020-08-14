#! /usr/bin/env python3
import argparse
import os

def main(args):
	names = []
	with open(args.list) as file:
		for index, line in enumerate(file):
			names.append(line.strip())
	print ("Loaded %s names from `%s`." % (len(names), args.list))
	commands = []
	for name in names:
		input_directory = args.inputs.rstrip("/")
		if args.subdirectory == 'y':
			sub_directory = name[:int(len(name)/2)]
			input_file = "%s/fasta/%s/%s.fasta" % (input_directory, sub_directory, name)
			output_path = "%s/hhr/%s" % (input_directory, sub_directory)
		else:
			input_file = "%s/fasta/%s.fasta" % (input_directory, name)
			output_path = "%s/hhr/" % (input_directory)
		if not os.path.exists(output_path):
			os.makedirs(output_path)
		output_file = "%s/%s.hhr" % (output_path, name)
		if os.path.isfile(output_file) and args.overwrite == 'n':
			print("Already available.")
			continue
		command = "%s -i %s -d %s -o %s" % (args.binary, input_file, args.database, output_file)
		os.system(command)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='This script performs HH-search on all matched files in the directory.')
	parser.add_argument('-l', '--list', help='Text file containing identifiers.', required=True)
	parser.add_argument('-i', '--inputs', help='Directory containing `fasta/X/Y.fasta` files', required=True)
	parser.add_argument('-b', '--binary', help='HH-search/HH-blits binary path', required=True)
	parser.add_argument('-d', '--database', help='HH-search database path', required=True)
	parser.add_argument('-s', '--subdirectory', help='Use subdirectory splitting', default='y')
	parser.add_argument('-w', '--overwrite', help='Overwrite existing results', default='n')
	args = parser.parse_args()
	main(args)