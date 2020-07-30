#! /usr/bin/env python3
import argparse
import os

def main(args):
	names = []
	with open(args.list) as file:
		for index, line in enumerate(file):
			names.append(line.strip())
	print ("Loaded %s names from `%s`." % (len(names), args.list))
	crossreference = {}
	with open(args.crossreference) as file:
		for index, line in enumerate(file):
			columns = line.split()
			core = columns[0]
			partner = columns[2]
			if core not in crossreference:
				crossreference[core] = []
			crossreference[core].append(partner)
	print ("Loaded cross reference from `%s`." % args.crossreference)
	targets = get_template_scores(args.target, args.minscore)
	print ("Loaded target scores from `%s`." % args.target)
	interactions = []
	for name in names:
		input_directory = args.inputs.rstrip("/")
		sub_directory = name[:int(len(name)/2)]			
		input_file = "%s/hhr/%s/%s.hhr" % (input_directory, sub_directory, name)
		templates = get_template_scores(input_file, args.minscore)
		minz = 0
		for t in targets:
			if t in crossreference:
				partners = crossreference[t]
				for p in partners:
					if p in templates:
						score = min(targets[t], templates[p])
						if score > minz:
							minz = score
		if minz > args.minscore:
			interactions.append((name, minz))
			print("Predicting: %s, min-Z: %s" % (name, minz))
	interactions.sort(key=lambda tup: tup[1], reverse=True)
	with open(args.output, 'w') as output_file:
		for i in interactions:
			output_file.write("%s %s\n" % (i[0], i[1]))

def get_template_scores(hhr_file, min_score):
	result = {}
	if os.path.isfile(hhr_file):
		with open(hhr_file) as file:
			for index, line in enumerate(file):
				if index > 8:
					if not line.strip():
						break
					template_id = line[4:10]
					template_score = float(line[57:63])
					if template_score > min_score:
						result[template_id] = template_score
	return result

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='This script identifies interactions by detecting matching HH-search results.')
	parser.add_argument('-t', '--target', help='HHR target file result', required=True)
	parser.add_argument('-c', '--crossreference', help='Cross Reference index file', required=True)
	parser.add_argument('-l', '--list', help='Text file containing identifiers.', required=True)
	parser.add_argument('-i', '--inputs', help='Directory containing `hhr/X/Y.hhr` files', required=True)
	parser.add_argument('-o', '--output', help='Output file containing minZ-scores`', required=True)
	parser.add_argument('-m', '--minscore', help='min-Z score threshold', default=10)
	args = parser.parse_args()
	main(args)