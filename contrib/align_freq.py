#!/usr/bin/python3

import re
import sys
import scipy


align_pair_re = re.compile(r"alignment for (\S+) and (\S+)")

# from http://stackoverflow.com/questions/10763697/histogram-without-plotting-function
def plot(data):
    """
    Histogram data to stdout
    """
    largest = max(data)
    scale = 50. / largest
    for i, datum in enumerate(data):
        bar = "*" * int(datum * scale)
        print ("%2d: %s (%d)" % (i, bar, datum))


def incr(lst, pos):
    if len(lst) <= pos:
        lst.extend([0 for i in range(pos - len(lst)  + 1)])
    lst[pos] += 1



if __name__ == "__main__":
    aligns = {}
    for line in open(sys.argv[1]):
        mobj = align_pair_re.search(line.strip())
        if mobj:
            target, query = mobj.groups()
            if not query in aligns: aligns[query] = []
            aligns[query].append(target)

    counts = []
    align_counts = []
    for k,v in aligns.items():
        print(k,len(v))
        counts.append(len(v))
        incr(align_counts, int(len(v)/10))

    print("Mean: ", scipy.mean(counts))
    print("StdDev: ", scipy.std(counts))
    print("Samples: ", len(counts))
    plot(align_counts)
