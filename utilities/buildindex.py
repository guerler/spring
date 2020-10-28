#! /usr/bin/env python3
import argparse
import os

def getSequences(fileName):
	sequences = dict()
	with open(fileName) as file:
		for line in file:
			if line.startswith(">"):
				name = line[1:7].upper()
				nextLine = next(file)
				sequences[name] = nextLine
	return sequences

def main(args):
	hhrEntries = set()
	with open(args.hhrlist) as file:
		for index, line in enumerate(file):
			hhrEntries.add(line[0:6].upper())
	print("Found %s hhr entries from `%s`." % (len(hhrEntries), args.hhrlist))

	dimers = set()
	with open(args.dimerlist) as file:
		for index, line in enumerate(file):
			dimers.add(line[0:4].upper())
	print("Found %s dimer entries from `%s`." % (len(dimers), args.dimerlist))

	hhrDimers = set()
	for entry in hhrEntries:
		if entry[0:4] in dimers:
			hhrDimers.add(entry)
	print("Found %s hhr entries in dimers." % len(hhrDimers))

	allEntries = dict()
	allCount = 0
	with open(args.fasta) as file:
		for index, line in enumerate(file):
			if line.startswith(">"):
				pdbId = line[1:5].upper()
				if pdbId not in allEntries:
					allEntries[pdbId] = set()
				allEntries[pdbId].add(line[1:7].upper())
				allCount = allCount + 1
	print("Found %s entries in complete fasta database." % allCount)

	crossReference = list()
	partnerList = set()
	for hhrDimer in hhrDimers:
		pdbId = hhrDimer[0:4]
		if pdbId not in allEntries:
			print("Warning: Missing entry %s" % pdbId)
			continue
		partners = allEntries[pdbId]
		if len(partners) > 0:
			for p in partners:
				crossReference.append([hhrDimer, p])
				if p not in hhrEntries:
					partnerList.add(p)
		else:
			crossReference.append([hhrDimer, hhrDimer])
	crossReference.sort(key=lambda x: (x[0], x[1]))
	print("Found %s index entries." % len(crossReference))
	print("Found %s additional binding partners for hhr entries." % len(partnerList))

	os.system("mkdir -p temp/")
	partnerListFile = "temp/partnerlist.txt"
	with open(partnerListFile, 'w') as output_file:
		for entry in partnerList:
			output_file.write("%s\n" % entry)

	partnerFasta = "temp/partner.fasta"
	os.system("./filterfasta.py -l %s -f %s -o %s" % (partnerListFile, args.fasta, partnerFasta))

	hhrFasta = "temp/hhr.fasta"
	os.system("./filterfasta.py -l %s -f %s -o %s" % (args.hhrlist, args.fasta, hhrFasta))
	os.system("makeblastdb -in %s -dbtype prot" % hhrFasta)

	hhrSequences = getSequences(hhrFasta)
	partnerSequences = getSequences(partnerFasta)
	print("Aligning %s on %s sequences..." % (len(partnerSequences.keys()), len(hhrSequences.keys())))

	partnerMatches = dict()
	partnerSequenceFile = "temp/query.fasta"
	partnerResultFile = "temp/query.out"
	partnerCount = 0
	for partnerEntry in partnerSequences.keys():
		maxScore = 0
		maxMatch = None
		partnerSequence = partnerSequences[partnerEntry]
		with open(partnerSequenceFile, 'w') as output_file:
			output_file.write(">%s\n" % partnerEntry)
			output_file.write("%s" % partnerSequence)
		os.system("psiblast -query %s -db %s -out %s" % (partnerSequenceFile, hhrFasta, partnerResultFile))
		try:
			with open(partnerResultFile) as file:
				for i in range(38):
					line = next(file)
				maxMatch = line.split()[0][0:6].upper()
		except:
			print("Warning: Skipping failed alignment [%s]." % partnerEntry)
			pass
		partnerMatches[partnerEntry] = maxMatch
		partnerCount = partnerCount + 1
		print("Matched %s: %s -> %s." % (partnerCount, partnerEntry, maxMatch))

	finalSet = set()
	for (core, partner) in crossReference:
		if partner in partnerList:
			if partner in partnerMatches:
				partner = partnerMatches[partner]
			else:
				partner = None
		if partner is not None:
			entry = "%s\t%s" % (core, partner)
			if entry not in finalSet:
				finalSet.add(entry)
		else:
			print("Warning: Skipping failed missing partner match [%s]." % partner)
	print("Found %s cross reference entries." % len(finalSet))

	with open(args.output, 'w') as output_file:
		for entry in sorted(finalSet):
			output_file.write("%s\n" % entry)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='List filtering.')
	parser.add_argument('-l', '--hhrlist', help='List of hhm template entries [PDB_CHAIN]', required=True)
	parser.add_argument('-d', '--dimerlist', help='List of all pdb multimers [PDB]', required=True)
	parser.add_argument('-f', '--fasta', help='Sequences of all pdb entries [PDB]', required=True)
	parser.add_argument('-o', '--output', help='Resulting list', required=True)
	args = parser.parse_args()
	main(args)