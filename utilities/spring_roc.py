#! /usr/bin/env python
import argparse
import math
import matplotlib.pyplot as plt
import os

def getId(rawId):
    elements = rawId.split("|")
    if len(elements) > 1:
        return elements[1]
    return rawId

def getReference(fileName, filterList=None, minScore=None, mappingDict=None, aCol=0, bCol=1, scoreCol=-1, separator=None):
    index = dict()
    failedToMap = 0
    with open(fileName) as fp:
        line = fp.readline()
        while line:
            ls = line.split(separator)
            a = getId(ls[aCol])
            b = getId(ls[bCol])
            skip = False
            if a == "-" or b == "-":
                skip = True
            aSplit = a.split("|")
            bSplit = b.split("|")
            if len(aSplit) > 1 or len(bSplit) > 1:
                print (aSplit)
                print (bSplit)
            if mappingDict is not None:
                if a in mappingDict and b in mappingDict:
                    a = mappingDict[a]
                    b = mappingDict[b]
                else:
                    failedToMap = failedToMap + 1
                    skip = True
            if filterList is not None:
                if a not in filterList and b not in filterList:
                    skip = True
            if not skip:
                if a > b:
                    name = "%s_%s" % (a, b)
                else:
                    name = "%s_%s" % (b, a)
                if name not in index:
                    if scoreCol >= 0 and len(ls) > scoreCol:
                        score = float(ls[scoreCol])
                        skip = False
                        if minScore is not None:
                            if minScore > score:
                                break
                        if not skip:
                            index[name] = score
                    else:
                        index[name] = 1.0
            line = fp.readline()
    return index

def getPercentage(rate, denominator):
    if denominator > 0:
        return 100.0 * rate / denominator
    return 0.0

def getxy(prediction, positive, negative):
    sorted_prediction = sorted(prediction.items(), key=lambda x: x[1], reverse=True)
    positive_total = len(positive)
    negative_total = len(negative)
    positive_denominator = float(positive_total)
    negative_denominator = float(negative_total)
    x = []
    y = []
    cnt = 0
    mcc = 0.0
    maxmcc = 0.0
    minscore = 0.0
    tp = 0
    fp = 0
    for (name, score) in sorted_prediction:
        cnt = cnt + 1
        found = False
        if name in positive:
            found = True
            tp = tp + 1
        if name in negative:
            found = True
            fp = fp + 1
        if found:
            x.append(getPercentage(fp, negative_denominator))
            y.append(getPercentage(tp, positive_denominator))
        fn = positive_total - tp
        tn = negative_total - fp
        denom = (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)
        if denom > 0.0:
            mcc = (tp*tn-fp*fn)/math.sqrt(denom)
            if mcc >= maxmcc:
                maxmcc = mcc
                minscore = score
    print("Top ranking prediction %s." % str(sorted_prediction[0]))
    print("Total count of prediction set: %s." % len(sorted_prediction))
    print("Total count of positive set: %s." % len(positive))
    print("Total count of negative set: %s." % len(negative))
    print("Matthews-Correlation-Coefficient: %s at Score >= %s." % (round(maxmcc, 2), minscore))
    return x, y

def main(args):
    # Check filter and mapping files
    mappingDict = None
    mappingName = "%s.mapping.txt" % args.input
    if os.path.isfile(mappingName):
        mappingDict = dict()
        with open(mappingName) as mappingFile:
            for line in mappingFile:
                ll = line.split()
                source = ll[0]
                target = ll[1]
                mappingDict[source] = target
        print("Total number in mapping list: %d." % len(mappingDict.keys()))

    filterList = None
    filterName = "%s.list.txt" % args.input
    if os.path.isfile(filterName):
        filterList = list()
        with open(filterName) as filterFile:
            for line in filterFile:
                id = getId(line.split()[0])
                if mappingDict is not None and id in mappingDict:
                    id = mappingDict[id]
                    filterList.append(id)
                else:
                    filterList.append(id)
        print("Total number in filter list: %d." % len(filterList))

    # Get data
    if args.biogrid:
        prediction = getReference("%s.txt" % args.input, scoreCol=2)
        positive = getReference(args.biogrid, filterList=filterList, aCol=23, bCol=26, separator="\t")
    else:
        prediction = getReference("%s.txt" % args.input, scoreCol=2, mappingDict=mappingDict)
        positive = getReference("%s.positive.txt" % args.input, filterList=filterList)
    negative = getReference("%s.negative.txt" % args.input)

    # Print plot
    x, y = getxy(prediction, positive, negative)
    plt.plot(x, y)
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