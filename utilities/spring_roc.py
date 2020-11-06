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

def getReference(fileName, filterList=None):
    index = dict()
    with open(fileName) as fp:
        line = fp.readline()
        while line:
            ls = line.split()
            a = getId(ls[0])
            b = getId(ls[1])
            skip = False
            if filterList is not None:
                if a not in filterList and b not in filterList:
                    skip = True
            if not skip and len(ls) >= 2:
                if a > b:
                    name = "%s_%s" % (a, b)
                else:
                    name = "%s_%s" % (b, a)
                if name not in index:
                    if len(ls) > 2:
                        index[name] = float(ls[2])
                    else:
                        index[name] = 1.0
            line = fp.readline()
    return index

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
            x.append(100.0 * fp / negative_denominator)
            y.append(100.0 * tp / positive_denominator)
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
    # Get filter list
    filterList = list()
    with open("%s.list" % args.input) as filterFile:
        for line in filterFile:
            filterList.append(getId(line))
    print(filterList)

    # Get data
    prediction = getReference("%s.txt" % args.input)
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
    args = parser.parse_args()
    main(args)