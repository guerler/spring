#! /usr/bin/env python
import argparse
import os

def main(args):
	files = os.listdir("%s" % args.inputs)
	for file in files:
		if file.endswith(".fasta"):
			name = os.path.splitext(file)[0]
			input_file = "%s/%s.fasta" % (args.inputs, name)
			output = "%s/%s.hhr" % (args.inputs, name)
			output_hhm = "%s/%s.hhm" % (args.inputs, name)
			output_a3m = "%s/%s.a3m" % (args.inputs, name)
			print("\nProcessing: %s\n" % (name))
			if os.path.isfile(output):
				print("Already available.\n")
				continue
			command = "%s -i %s -d %s -ohhm %s -oa3m %s &" % (args.binary, input_file, args.database, output_hhm, output_a3m)
			print ("Executing: %s\n" % command)
			raw_input('Press ENTER to continue or CTRL-C to exit.')
			os.system(command)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='This script performs HH-search on all matched files in the directory.')
	parser.add_argument('-i', '--inputs', help='Directory with files *.fasta', required=True)
	parser.add_argument('-b', '--binary', help='HH-search/HH-blits binary path', required=True)
	parser.add_argument('-d', '--database', help='HH-search database path', required=True)
	args = parser.parse_args()
	main(args)