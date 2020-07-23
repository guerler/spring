#! /usr/bin/env python
import argparse
import os

def main(args):
	names = []
	with open(args.list) as file:
		for index, line in enumerate(file):
			names.append(line.strip())
	print ("Loaded %s names from `%s`." % (len(names), args.list))
	for name in names:
		print("\nProcessing: %s" % name)
		print(30 * "-")
		input_directory = args.inputs.rstrip("/")
		sub_directory = name[:len(name)/2]			
		input_file = "%s/fasta/%s/%s.fasta" % (input_directory, sub_directory, name)
		output_path = "%s/hhr/%s" % (input_directory, sub_directory)
		if not os.path.exists(output_path):
			os.makedirs(output_path)
		output_file = "%s/%s.hhr" % (output_path, name)
		if os.path.isfile(output_file):
			print("Already available.")
			continue
		background_flag = "" if args.fg else "&";
		command = "%s -i %s -d %s -ohhr %s %s" % (args.binary, input_file, output_file, args.database, background_flag)
		print ("Executing: %s" % command)
		if not args.enter:
			raw_input('Press ENTER to continue or CTRL-C to exit.')
		os.system(command)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='This script performs HH-search on all matched files in the directory.')
	parser.add_argument('-l', '--list', help='Text file containing identifiers.', required=True)
	parser.add_argument('-i', '--inputs', help='Directory containing `fasta/X/Y.fasta` files', required=True)
	parser.add_argument('-b', '--binary', help='HH-search/HH-blits binary path', required=True)
	parser.add_argument('-d', '--database', help='HH-search database path', required=True)
	parser.add_argument('--enter', help='Do not ask user to confirm each command.', action='store_true')
	parser.add_argument('--fg', help='Disable background submission.', action='store_true')
	args = parser.parse_args()
	main(args)