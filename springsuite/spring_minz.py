#! /usr/bin/env python3
import argparse
import os

def main(args):
	inputs = list()
	with open(args.inputlist) as file:
		for index, line in enumerate(file):
			name = line.strip()
			inputs.append(name)
	print ("Loaded %s input names from `%s`." % (len(inputs), args.inputlist))
	targets = list()
	duplicates = 0
	with open(args.targetlist) as file:
		for index, line in enumerate(file):
			name = line.strip()
			targets.append(name)
			if name in inputs:
				duplicates = duplicates + 1
	print ("Loaded %s target names from `%s`." % (len(targets), args.targetlist))
	crossReference = dict()
	with open(args.crossreference) as file:
		for index, line in enumerate(file):
			columns = line.split()
			core = columns[0]
			partner = columns[-1]
			if core not in crossReference:
				crossReference[core] = []
			crossReference[core].append(partner)
	print ("Loaded cross reference from `%s`." % args.crossreference)
	interactions = dict()
	for targetName in targets:
		targetDirectory = args.targetpath.rstrip("/")
		targetFile = "%s/%s" % (targetDirectory, targetName)
		matchScores(targetFile=targetFile,
					targetName=targetName,
					inputs=sorted(inputs),
					inputPath=args.inputpath,
					crossReference=crossReference,
					minScore=args.minscore,
					idLength=args.idx,
					interactions=interactions)
	if duplicates != len(targets):
		for inputName in inputs:
			inputDirectory = args.inputpath.rstrip("/")
			inputFile = "%s/%s" % (inputDirectory, inputName)
			matchScores(targetFile=inputFile,
						targetName=inputName,
						inputs=targets,
						inputPath=args.targetpath,
						crossReference=crossReference,
						minScore=args.minscore,
						idLength=args.idx,
						interactions=interactions)
	interactions = sorted(interactions.values(), key=lambda item: item["minZ"], reverse=True)
	with open(args.output, 'w') as output_file:
		for entry in interactions:
			output_file.write("%s\t%s\t%s\t%s\n" % (entry["targetName"], entry["inputName"], entry["minZ"], entry["minInfo"]))

def matchScores(targetFile, targetName, inputs, inputPath, crossReference, minScore, idLength, interactions):
	targetTop, targetHits = getTemplateScores(targetFile, minScore, idLength)
	if not targetHits:
		print("No targets found `%s`" % targetFile)
	else:
		print ("Loaded target scores from `%s`." % targetFile)
		for inputName in inputs:
			inputDirectory = inputPath.rstrip("/")
			inputFile = "%s/%s" % (inputDirectory, inputName)
			inputTop, inputHits = getTemplateScores(inputFile, minScore, idLength)
			minZ = 0
			minInfo = ""
			for t in targetHits:
				if t in crossReference:
					partners = crossReference[t]
					for p in partners:
						if p in inputHits:
							score = min(targetHits[t], inputHits[p])
							if score > minZ:
								minZ = score
								minInfo = "%s\t%s\t%s\t%s" % (targetTop, inputTop, t, p)
			if minZ > minScore:
				if targetName > inputName:
					interactionKey = "%s_%s" % (targetName, inputName)
				else:
					interactionKey = "%s_%s" % (inputName, targetName)
				if interactionKey in interactions:
					if interactions[interactionKey]["minZ"] >= minZ:
						continue
				interactions[interactionKey] = dict(targetName=targetName, inputName=inputName, minZ=minZ, minInfo=minInfo)
				print("Predicting: %s, min-Z: %s, templates: %s" % (inputName, minZ, minInfo))
	return interactions

def getTemplateScores(hhrFile, minScore, identifierLength):
	result = dict()
	topTemplate = None
	identifierLength = identifierLength + 4
	if os.path.isfile(hhrFile):
		with open(hhrFile) as file:
			for index, line in enumerate(file):
				if index > 8:
					if not line.strip():
						break
					templateId = line[4:identifierLength]
					templateScore = float(line[57:63])
					if templateScore > minScore:
						if topTemplate is None:
							topTemplate = templateId
						result[templateId] = templateScore
	return topTemplate, result

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='This script identifies interactions by detecting matching HH-search results.')
	parser.add_argument('-il', '--inputlist', help='Text file containing identifiers.', required=True)
	parser.add_argument('-ip', '--inputpath', help='Directory containing `hhr` files', required=True)
	parser.add_argument('-tl', '--targetlist', help='Text file containing identifiers.', required=True)
	parser.add_argument('-tp', '--targetpath', help='Directory containing `hhr` files', required=True)
	parser.add_argument('-c', '--crossreference', help='Cross Reference index file', required=True)
	parser.add_argument('-x', '--idx', help='Length of identifier', type=int, default=6)
	parser.add_argument('-o', '--output', help='Output file containing min-Z scores', required=True)
	parser.add_argument('-m', '--minscore', help='min-Z score threshold', type=int, default=10)
	args = parser.parse_args()
	main(args)