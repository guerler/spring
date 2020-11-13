#! /usr/bin/env python3
import argparse
import os


def main(args):
    logFile = open(args.log, 'a+')
    targets = list()
    targetPath = args.targetpath.rstrip("/")
    with open(args.targetlist) as file:
        for line in file:
            name = line.strip()
            targets.append(name+ ".txt")
    print("Loaded %s target names from `%s`." % (len(targets),
          args.targetlist))
    if args.inputlist:
        inputs = list()
        inputPath = args.inputpath.rstrip("/")
        with open(args.inputlist) as file:
            for line in file:
                name = line.strip()
                inputs.append(name )
        print("Loaded %s input names from `%s`." % (len(inputs),
              args.inputlist))
    else:
        inputs = targets
        inputPath = targetPath
    crossReference = dict()
    with open(args.crossreference) as file:
        for line in file:
            columns = line.split()
            core = columns[0]
            partner = columns[-1]
            if core not in crossReference:
                crossReference[core] = []
            crossReference[core].append(partner)
    print("Loaded cross reference from `%s`." % args.crossreference)
    interactions = dict()
    for targetName in targets:
        targetFile = "%s/%s" % (targetPath, targetName)
        matchScores(targetFile=targetFile,
                    targetName=targetName,
                    inputs=inputs,
                    inputPath=inputPath,
                    crossReference=crossReference,
                    idLength=args.idlength,
                    minScore=args.minscore,
                    logFile=logFile,
                    interactions=interactions)
    if args.inputlist:
        for inputName in inputs:
            inputDirectory = inputPath
            inputFile = "%s/%s" % (inputDirectory, inputName)
            matchScores(targetFile=inputFile,
                        targetName=inputName,
                        inputs=targets,
                        inputPath=targetPath,
                        crossReference=crossReference,
                        minScore=args.minscore,
                        idLength=args.idlength,
                        logFile=logFile,
                        interactions=interactions)
    interactions = sorted(interactions.values(), key=lambda item: item["minZ"],
                          reverse=True)
    with open(args.output, 'w') as output_file:
        for entry in interactions:
            output_file.write("%s\t%s\t%s\t%s\n" % (entry["targetName"],
                              entry["inputName"], entry["minZ"],
                              entry["minInfo"]))
    logFile.close()


def matchScores(targetFile, targetName, inputs, inputPath, crossReference,
                minScore, idLength, logFile, interactions):
    targetTop, targetHits = getTemplateScores(targetFile, minScore, idLength)
    if not targetHits:
        print("No targets found `%s`" % targetFile)
    else:
        print("Loaded target scores from `%s`." % targetFile)
        for inputName in inputs:
            inputFile = "%s/%s" % (inputPath, inputName)
            inputTop, inputHits = getTemplateScores(inputFile,
                                                    minScore, idLength)
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
                                minInfo = "%s\t%s\t%s\t%s" % (targetTop,
                                                              inputTop, t, p)
            if minZ > minScore:
                if targetName > inputName:
                    interactionKey = "%s_%s" % (targetName, inputName)
                else:
                    interactionKey = "%s_%s" % (inputName, targetName)
                if interactionKey in interactions:
                    if interactions[interactionKey]["minZ"] >= minZ:
                        continue
                interactions[interactionKey] = dict(targetName=targetName,
                                                    inputName=inputName,
                                                    minZ=minZ, minInfo=minInfo)
                logFile.write("Interaction between %s and %s [min-Z: %s].\n" % 
                              (targetName, inputName, minZ))


def getTemplateScores(hhrFile, minScore, idLength):
    result = dict()
    topTemplate = None
    idLength = idLength + 4
    if os.path.isfile(hhrFile):
        with open(hhrFile) as file:
            for index, line in enumerate(file):
                if index > 8:
                    if not line.strip():
                        break
                    templateId = line[4:idLength]
                    templateScore = float(line[57:63])
                    if templateScore > minScore:
                        if topTemplate is None:
                            topTemplate = templateId
                        result[templateId] = templateScore
    return topTemplate, result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script identifies interactions by detecting matching HH-search results.')
    parser.add_argument('-tl', '--targetlist', help='Text file containing identifiers.', required=True)
    parser.add_argument('-tp', '--targetpath', help='Directory containing `hhr` files', required=True)
    parser.add_argument('-il', '--inputlist', help='Text file containing identifiers.', required=False)
    parser.add_argument('-ip', '--inputpath', help='Directory containing `hhr` files', required=False)
    parser.add_argument('-c', '--crossreference', help='Cross Reference index file', required=True)
    parser.add_argument('-o', '--output', help='Output file containing min-Z scores', required=True)
    parser.add_argument('-l', '--log', help='Log file', required=True)
    parser.add_argument('-m', '--minscore', help='min-Z score threshold', type=int, default=10)
    parser.add_argument('-idx', '--idlength', help='Length of identifier in reference', type=int, default=6)
    args = parser.parse_args()
    main(args)
