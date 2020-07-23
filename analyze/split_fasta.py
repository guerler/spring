#! /usr/bin/env python
import argparse
import os

def main(args):
	names = []
	with open(args.inputs) as file:
		current_name = None
		current_content = ''
		for index, line in enumerate(file):
			if line.startswith('>'):
				if current_name and len(current_name) > 2:
					output_subdirectory = current_name[:len(current_name)/2]
					output_path = "%s/%s" % (args.output.rstrip("/"), output_subdirectory)
					if not os.path.exists(output_path):
						os.makedirs(output_path)
					output_name = "%s/%s.fasta" % (output_path, current_name)
					names.append(current_name)
					with open(output_name, 'w') as output_file:
						output_file.write(">%s\n" % current_name)
						output_file.write(current_content)
					current_content = ''
				current_name = line[1:-1].split(' ', 1)[0]
			elif current_name is not None:
				current_content = current_content + line
	output_index_file = "output.txt"
	print ("Created %s fasta files. Names written to '%s'." % (len(names), output_index_file))
	with open(output_index_file, 'w') as output_index:
		output_index.write("\n".join(names))

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='This script splits a file containing multiple fasta sequences into individual files.')
	parser.add_argument('-i', '--inputs', help='File containing fasta sequences', required=True)
	parser.add_argument('-o', '--output', help='Path to output directory', required=True)
	args = parser.parse_args()
	main(args)