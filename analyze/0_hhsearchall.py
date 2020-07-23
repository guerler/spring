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
			print("\nProcessing: %s" % name)
			print(30 * "-")
			if os.path.isfile(output):
				print("Already available.")
				continue
			background_flag = "" if args.fg else "&";
			command = "%s -i %s -d %s -ohhm %s -oa3m %s %s" % (args.binary, input_file, args.database, output_hhm, output_a3m, background_flag)
			print ("Executing: %s" % command)
			if not args.enter:
				raw_input('Press ENTER to continue or CTRL-C to exit.')
			os.system(command)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='This script performs HH-search on all matched files in the directory.')
	parser.add_argument('-i', '--inputs', help='Directory with files *.fasta', required=True)
	parser.add_argument('-b', '--binary', help='HH-search/HH-blits binary path', required=True)
	parser.add_argument('-d', '--database', help='HH-search database path', required=True)
	parser.add_argument('--enter', help='Do not ask user to confirm each command.', action='store_true')
	parser.add_argument('--fg', help='Disable background submission.', action='store_true')
	args = parser.parse_args()
	main(args)