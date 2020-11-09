#! /usr/bin/env python
import argparse
import math
import matplotlib.pyplot as plt
import os

def getIds(rawIds):
    return rawIds.split("|")

def getCenterId(rawId):
    elements = rawId.split("|")
    if len(elements) > 1:
        return elements[1]
    return rawId

def getReference(fileName, filterList=None, minScore=None, aCol=0, bCol=1, scoreCol=-1, separator=None, skipStartsWith="SWISS-PROT"):
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
            for a in aList:
                for b in bList:
                    skip = False
                    if a == "-" or b == "-" or a.startswith(skipStartsWith):
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
            count = count + 1
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
    count = 0
    mcc = 0.0
    maxcount = 0
    maxmcc = 0.0
    maxprecision = 0.0
    minscore = 0.0
    tp = 0
    fp = 0
    previous = 0.0
    for (name, score) in sorted_prediction:
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
                maxcount = count
                maxprecision = 0.0
                if tp > 0 or fp > 0:
                    maxprecision = tp / (tp + fp)
        if count % 10000 == 0:
            print ("%s (precision=%5.3f)" % (count, maxprecision))
        count = count + 1

    print("Top ranking prediction %s." % str(sorted_prediction[0]))
    print("Total count of prediction set: %s (Precision=%1.2f)." % (maxcount, maxprecision))
    print("Total count of positive set: %s." % len(positive))
    print("Total count of negative set: %s." % len(negative))
    print("Matthews-Correlation-Coefficient: %s at Score >= %s." % (round(maxmcc, 2), minscore))
    return x, y

def main(args):
    # filter list
    filterList = None
    filterName = "%s.list.txt" % args.input
    if os.path.isfile(filterName):
        filterList = set()
        with open(filterName) as filterFile:
            for line in filterFile:
                id = getCenterId(line.split()[0])
                filterList.add(id)
        print("Total number in filter list: %d." % len(filterList))

    # process biogrid database
    if args.biogrid:
        positive = getReference(args.biogrid, aCol=23, bCol=26, separator="\t", filterList=filterList)
        prediction = getReference("%s.txt" % args.input, scoreCol=2)
    else:
        positive = getReference("%s.positive.txt" % args.input, filterList=filterList)
        prediction = getReference("%s.txt" % args.input, scoreCol=2)

    # estimate background noise
    negative = getReference("%s.negative.txt" % args.input)

    # print plot
    print ("Producing plot data...")
    print("Total count in prediction file: %d." % len(prediction))
    print("Total count in positive file: %d." % len(positive))
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