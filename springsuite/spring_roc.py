#! /usr/bin/env python
import argparse
import math
import random
from datetime import datetime

from matplotlib import pyplot as plt


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
                 skipFirstLine=False, filterValues=list()):
    index = dict()
    count = 0
    with open(fileName) as fp:
        line = fp.readline()
        if skipFirstLine:
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
                    if a == "-" or b == "-":
                        skip = True
                    if filterA is not None:
                        if a not in filterA and b not in filterA:
                            skip = True
                    if filterB is not None:
                        if a not in filterB and b not in filterB:
                            skip = True
                    for f in filterValues:
                        if len(ls) > f[0]:
                            columnEntry = ls[f[0]].lower()
                            searchEntry = f[1].lower()
                            if columnEntry.find(searchEntry) == -1:
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
    x = list([0])
    y = list([0])
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
            yValue = getPercentage(tp, tp + fn)
            xValue = getPercentage(fp, fp + tn)
            y.append(yValue)
            x.append(xValue)
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
    filterSets = getFilter(args.input)
    filterKeys = list(filterSets.keys())
    filterA = filterSets[filterKeys[0]]
    if len(filterKeys) > 1:
        filterB = filterSets[filterKeys[1]]
    else:
        filterB = filterA

    # identify biogrid filter options
    filterValues = []
    if args.method:
        filterValues.append([11, args.method])
    if args.experiment:
        filterValues.append([12, args.experiment])
    if args.throughput:
        filterValues.append([17, args.throughput])

    # process biogrid database
    print("Loading positive set from BioGRID file...")
    positive, positiveCount = getReference(args.biogrid, aCol=23, bCol=26,
                                           separator="\t", filterA=filterA,
                                           filterB=filterB, skipFirstLine=True,
                                           filterValues=filterValues)

    # rescan biogrid database to identify set of putative interactions
    if filterValues:
        print("Filtered entries by (column, value): %s" % filterValues)
        print("Loading putative set from BioGRID file...")
        putative, putativeCount = getReference(args.biogrid, aCol=23, bCol=26,
                                               separator="\t", filterA=filterA,
                                               filterB=filterB,
                                               skipFirstLine=True)
        print("Found %s." % putativeCount)
    else:
        putative = positive

    # process prediction file
    print("Loading prediction file...")
    prediction, _ = getReference(args.input, scoreCol=2)

    # check if localization file has been provided otherwise generate random background
    if args.localization:
        # get subcellular localization from UniProt export
        localA = dict()
        localB = dict()
        with open(args.localization) as locFile:
            for line in locFile:
                searchKey = "SUBCELLULAR LOCATION"
                searchPos = line.find(searchKey)
                if searchPos != -1:
                    uniId = line.split()[0]
                    locStart = searchPos + len(searchKey) + 1
                    locId = line[locStart:].split()[0]
                    if locId in ["Nucleus", "Membrane", "Cytoplasm"]:
                        if uniId in filterA:
                            localA[uniId] = locId
                        if uniId in filterB:
                            localB[uniId] = locId

        localAKeys = list(localA.keys())
        localBKeys = list(localB.keys())
        print("Found %d, %d subcellular localizations." % (len(localAKeys), len(localBKeys)))

        # estimate background noise
        print("Estimating background noise...")
        negative = set()
        negativeRequired = positiveCount
        for nameA in localAKeys:
            for nameB in localBKeys:
                if localA[nameA] != localB[nameB]:
                    key = getKey(nameA, nameB)
                    if key not in putative and key not in negative and negativeRequired > 0:
                        negative.add(key)
                        negativeRequired = negativeRequired - 1
            if negativeRequired <= 0:
                break
    else:
        # estimate background noise
        print("Estimating background noise...")
        negative = set()
        filterAList = list(filterA)
        filterBList = list(filterB)
        negativeCount = positiveCount
        negativeRequired = negativeCount
        random.seed(datetime.now())
        while negativeRequired > 0:
            nameA = random.choice(filterAList)
            nameB = random.choice(filterBList)
            key = getKey(nameA, nameB)
            if key not in putative and key not in negative:
                negative.add(key)
                negativeRequired = negativeRequired - 1

    # create plot
    print("Producing plot data...")
    print("Total count in prediction file: %d." % len(prediction))
    print("Total count in positive file: %d." % len(positive))
    plt.ylabel('True Positive Rate (%)')
    plt.xlabel('False Positive Rate (%)')
    title = " vs. ".join(filterSets)
    plt.suptitle(title)
    if filterValues:
        filterAttributes = list(map(lambda x: x[1], filterValues))
        plt.title("BioGRID filters: %s" % filterAttributes, fontsize=10)
    x, y, xMax = getXY(prediction, positive, positiveCount, negative)
    plt.plot(x, y)
    plt.plot([0, xMax], [0, xMax])
    plt.savefig(args.output, format="png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create ROC plot.')
    parser.add_argument('-i', '--input', help='Input prediction file.', required=True)
    parser.add_argument('-b', '--biogrid', help='BioGRID interaction database file', required=True)
    parser.add_argument('-l', '--localization', help='UniProt export table with subcellular localizations', required=False)
    parser.add_argument('-e', '--experiment', help='Type (physical/genetic)', required=False)
    parser.add_argument('-t', '--throughput', help='Throughput (low/high)', required=False)
    parser.add_argument('-m', '--method', help='Method e.g. Two-hybrid', required=False)
    parser.add_argument('-o', '--output', help='Output (png)', required=True)
    args = parser.parse_args()
    main(args)
