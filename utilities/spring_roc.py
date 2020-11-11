#! /usr/bin/env python
import argparse
import math
from matplotlib import pyplot as plt
import random
from datetime import datetime


def getIds(rawIds):
    return rawIds.split("|")


def getCenterId(rawId):
    elements = rawId.split("|")
    if len(elements) > 1:
        return elements[1]
    return rawId


def getOrganism(rawId):
    elements = rawId.split("_")
    return elements[-1]


def getKey(a, b):
    if a > b:
        name = "%s_%s" % (a, b)
    else:
        name = "%s_%s" % (b, a)
    return name


def getPercentage(rate, denominator):
    if denominator > 0:
        return 100.0 * rate / denominator
    return 0.0


def getFilter(filterName):
    print("Loading target organism(s)...")
    filterSets = dict()
    with open(filterName) as filterFile:
        for line in filterFile:
            columns = line.split()
            for colIndex in [0, 1]:
                if colIndex >= len(columns):
                    break
                colEntry = columns[colIndex]
                id = getCenterId(colEntry)
                organism = getOrganism(colEntry)
                if organism not in filterSets:
                    filterSets[organism] = set()
                filterSets[organism].add(id)
    print("Organism(s) in set: %s." % filterSets.keys())
    return filterSets


def getReference(fileName, filterA=None, filterB=None, minScore=None, aCol=0,
                 bCol=1, scoreCol=-1, separator=None,
                 skipStartsWith="SWISS-PROT"):
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
                bList = [getCenterId(ls[bCol])]
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


def getXY(prediction, positive, positiveCount, negative):
    sortedPrediction = sorted(prediction.items(), key=lambda x: x[1],
                              reverse=True)
    positiveTotal = positiveCount
    negativeTotal = len(negative)
    positiveDenominator = float(positiveTotal)
    negativeDenominator = float(negativeTotal)
    x = list()
    y = list()
    xMax = 0
    topCount = 0
    topMCC = 0.0
    topPrecision = 0.0
    topScore = 0.0
    tp = 0
    fp = 0
    count = 0
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
            if mcc >= topMCC:
                topMCC = mcc
                topScore = score
                topCount = count
                topPrecision = precision
        if found:
            yValue = getPercentage(tp, positiveDenominator)
            xValue = getPercentage(fp, negativeDenominator)
            x.append(xValue)
            y.append(yValue)
            xMax = max(xValue, xMax)
        count = count + 1
    print("Top ranking prediction %s." % str(sortedPrediction[0]))
    print("Total count of prediction set: %s (precision=%1.2f)." %
          (topCount, topPrecision))
    print("Total count of positive set: %s." % len(positive))
    print("Total count of negative set: %s." % len(negative))
    print("Matthews-Correlation-Coefficient: %s at Score >= %s." %
          (round(topMCC, 2), topScore))
    return x, y, xMax


def main(args):
    # load source files
    filterSets = getFilter("%s.list" % args.input)
    filterKeys = list(filterSets.keys())
    filterA = filterSets[filterKeys[0]]
    if len(filterKeys) > 1:
        filterB = filterSets[filterKeys[1]]
    else:
        filterB = filterA

    # process biogrid database
    print("Loading postive set from BioGRID file...")
    positive, positiveCount = getReference(args.biogrid, aCol=23, bCol=26,
                                           separator="\t", filterA=filterA,
                                           filterB=filterB)

    # process prediction file
    print("Loading prediction file...")
    prediction, _ = getReference(args.input, scoreCol=2)

    # estimate background noise
    print("Estimating background noise...")
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

    # create plot
    print("Producing plot data...")
    print("Total count in prediction file: %d." % len(prediction))
    print("Total count in positive file: %d." % len(positive))
    x, y, xMax = getXY(prediction, positive, positiveCount, negative)
    plt.plot(x, y)
    plt.plot([0, xMax], [0, xMax])
    plt.xlabel('False Positive Rate (%)')
    plt.ylabel('True Positive Rate (%)')
    plt.title('Prediction Validation')
    plt.savefig("out.pdf", formatstr="pdf")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create ROC plot.')
    parser.add_argument('-i', '--input', help='Input prediction file.',
                        required=True)
    parser.add_argument('-b', '--biogrid', help='BioGRID interaction ' +
                        'database file', required=True)
    args = parser.parse_args()
    main(args)
