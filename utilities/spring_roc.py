#! /usr/bin/env python
import argparse
import math
import matplotlib.pyplot as plt
import os
import random
from datetime import datetime

def getIds(rawIds):
    return rawIds.split("|")

def getCenterId(rawId):
    elements = rawId.split("|")
    if len(elements) > 1:
        return elements[1]
    return rawId

def getKey(a, b):
    if a > b:
        name = "%s_%s" % (a, b)
    else:
        name = "%s_%s" % (b, a)
    return name

def getReference(fileName, filterA=None, filterB=None, minScore=None, aCol=0, bCol=1, scoreCol=-1, separator=None, skipStartsWith="SWISS-PROT"):
    index = dict()
    count = 0
    with open(fileName) as fp:
        line = fp.readline()
        while line:
            ls = line.split(separator)
            if separator is not None:
                aList = getIds(ls[aCol])
                bList = getIds(ls[bCol])
            else:
                aList = [getCenterId(ls[aCol])]
                bList= [getCenterId(ls[bCol])]
            validEntry = False
            for a in aList:
                for b in bList:
                    skip = False
                    if a == "-" or b == "-" or a.startswith(skipStartsWith):
                        skip = True
                    if filterA is not None:
                        if a not in filterA and b not in filterA:
                            skip = True
                    if filterB is not None:
                        if a not in filterB and b not in filterB:
                            skip = True
                    if not skip:
                        name = getKey(a, b)
                        if name not in index:
                            validEntry = True
                            if scoreCol >= 0 and len(ls) > scoreCol:
                                score = float(ls[scoreCol])
                                skip = False
                                if minScore is not None:
                                    if minScore > score:
                                        return index, count
                                if not skip:
                                    index[name] = score
                            else:
                                index[name] = 1.0
            if validEntry:
                count = count + 1
            line = fp.readline()
    return index, count

def getPercentage(rate, denominator):
    if denominator > 0:
        return 100.0 * rate / denominator
    return 0.0

def getXY(prediction, positive, positiveCount, negative):
    sortedPrediction = sorted(prediction.items(), key=lambda x: x[1], reverse=True)
    positiveTotal = positiveCount
    negativeTotal = len(negative)
    positiveDenominator = float(positiveTotal)
    negativeDenominator = float(negativeTotal)
    x = []
    y = []
    xMax = 0
    count = 0
    mcc = 0.0
    maxcount = 0
    maxmcc = 0.0
    maxprecision = 0.0
    maxscore = 0.0
    tp = 0
    fp = 0
    for (name, score) in sortedPrediction:
        found = False
        if name in positive:
            found = True
            tp = tp + 1
        if name in negative:
            found = True
            fp = fp + 1
        precision = 0.0
        if tp > 0 or fp > 0:
            precision = tp / (tp + fp)
        fn = positiveTotal - tp
        tn = negativeTotal - fp
        denom = (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)
        if denom > 0.0:
            mcc = (tp*tn-fp*fn)/math.sqrt(denom)
            if mcc >= maxmcc:
                maxmcc = mcc
                maxscore = score
                maxcount = count
                maxprecision = precision
        if found:
            yValue = getPercentage(tp, positiveDenominator)
            xValue = getPercentage(fp, negativeDenominator)
            x.append(xValue)
            y.append(yValue)
            xMax = xValue if xValue > xMax else xMax
        if count % 10000 == 0:
            print ("%s (precision=%5.3f)" % (count, maxprecision))
        count = count + 1

    print("Top ranking prediction %s." % str(sortedPrediction[0]))
    print("Total count of prediction set: %s (precision=%1.2f)." % (maxcount, maxprecision))
    print("Total count of positive set: %s." % len(positive))
    print("Total count of negative set: %s." % len(negative))
    print("Matthews-Correlation-Coefficient: %s at Score >= %s." % (round(maxmcc, 2), maxscore))
    return x, y, xMax + 1

def getFilter(filterName):
    filterSet = set()
    if os.path.isfile(filterName):
        with open(filterName) as filterFile:
            for line in filterFile:
                id = getCenterId(line.split()[0])
                filterSet.add(id)
    print("Total number in filter list: %d." % len(filterSet))
    return filterSet

def main(args):
    inputBase = os.path.basename(args.input)
    inputPath = os.path.dirname(os.path.abspath(args.input))
    inputName = os.path.splitext(inputBase)[0]
    splitTargets = inputName.split("_")
    print("Target names: %s." % splitTargets)
    if len(splitTargets) > 1:
        filterA = getFilter("%s/%s.list.txt" % (inputPath, splitTargets[0]))
        filterB = getFilter("%s/%s.list.txt" % (inputPath, splitTargets[1]))
    else:
        filterA = filterB = getFilter("%s/%s.list.txt" % (inputPath, splitTargets[0]))

    # process biogrid database
    if args.biogrid:
        positive, positiveCount = getReference(args.biogrid, aCol=23, bCol=26, separator="\t", filterA=filterA, filterB=filterB)
        prediction, _ = getReference("%s.txt" % args.input, scoreCol=2)
    else:
        positive, positiveCount = getReference("%s.positive.txt" % args.input, filterA=filterA, filterB=filterB)
        prediction, _ = getReference("%s.txt" % args.input, scoreCol=2)

    # estimate background noise
    negative = set()
    filterAList = list(filterA)
    filterBList = list(filterB)
    negativeCount = int(len(positive.keys()))
    random.seed(datetime.now())
    while negativeCount > 0:
        nameA = random.choice(filterAList)
        nameB = random.choice(filterBList)
        key = getKey(nameA, nameB)
        if key not in positive and key not in negative:
            negative.add(key)
            negativeCount = negativeCount - 1

    # print plot
    print ("Producing plot data...")
    print("Total count in prediction file: %d." % len(prediction))
    print("Total count in positive file: %d." % len(positive))
    x, y, xMax = getXY(prediction, positive, positiveCount, negative)
    plt.plot(x, y)
    plt.plot(range(int(xMax)), range(int(xMax)))
    plt.xlabel('False Positive Rate (%)')
    plt.ylabel('True Positive Rate (%)')
    plt.title('Prediction Results')
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create ROC plot.')
    parser.add_argument('-i', '--input', help='Input template entries [PDB_CHAIN]', required=True)
    parser.add_argument('-b', '--biogrid', help='Read BioGrid interaction database file', required=False)
    args = parser.parse_args()
    main(args)