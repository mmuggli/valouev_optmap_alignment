import re
import sys
import scipy
import scipy.stats
align_pair_re = re.compile(r"alignment for (\S+) and (\S+)")
alignment_line_re = re.compile(r"\[(.*)\]->\[(.*)\] s: (.*)")
frag_re = re.compile(r"(\d+):(\d+(\.\d+)?)")


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


def sumfrags(frags):
    mysum = 0.0
    for frag in frags:
        mysum += float(frag[1])
    return mysum

if __name__ == "__main__":
    matches = []
    num_frags = []
    nontriv_matches = []
    alignswith3 = []
    dofs = []
    chisquaredsums = []
    STDDEV = 2458
    for line in open(sys.argv[1]):
        if align_pair_re.search(line):
            matches.append(0)
            nontriv_matches.append(0)
            alignswith3.append(0)
            dofs.append(0)
            chisquaredsums.append(0.0)
            
        mobj = alignment_line_re.search(line)
        if mobj:
            frags_l = frag_re.findall(mobj.group(1))
            frags_r = frag_re.findall(mobj.group(2))

            # chi squared stuff
            #if len(frags_l) == 1 and len(frags_r) == 1:
            if True:
                diff = (sumfrags(frags_l) - sumfrags(frags_r)) * 1000
                chisquaredsums[-1] += (abs(diff) / STDDEV) ** 2
                dofs[-1] += 1
                    
            incr(num_frags, len(frags_l))
            incr(num_frags,len(frags_r))
            if len(frags_l) > 1:
                nontriv_matches[-1] += 1
            if len(frags_l) > 2:
                alignswith3[-1] = +1

        matches[-1] += 1

    print("--- Number of matched fragments or combined fragments per alignment ---")
    print("Mean: ", scipy.mean(matches))
    print("StdDev: ", scipy.std(matches))
    print("Samples: ", len(matches))
    print("Max: ", max(matches))
    print("Min: ", min(matches))

    print("--- Number of non-trivial matched combined fragments on left side per alignment ---")
    print("Mean: ", scipy.mean(nontriv_matches))
    print("StdDev: ", scipy.std(nontriv_matches))
    print("Samples: ", len(nontriv_matches))
    print("Max: ", max(nontriv_matches))
    print("Min: ", min(nontriv_matches))


    print("--- Number frags per possibly combined fragment ---")
    print("Mean: ", scipy.mean(num_frags))
    print("StdDev: ", scipy.std(num_frags))
    print("Samples: ", len(num_frags))
    print("Max: ", max(num_frags))
    print("Min: ", min(num_frags))
    plot(num_frags)


    print("--- Number of three frag matches per alignment ---")
    print("Mean: ", scipy.mean(alignswith3))
    print("StdDev: ", scipy.std(alignswith3))
    print("Samples: ", len(alignswith3))
    print("Max: ", max(alignswith3))
    print("Min: ", min(alignswith3))
    print("nonzero count: ", sum([1 for i in alignswith3 if i > 0]))


    chisquareds = [scipy.stats.chi2.cdf(chisquaredsum, dof) for chisquaredsum, dof in zip(chisquaredsums, dofs)]
    print("--- chi squared cdfs ---")
    print("Mean: ", scipy.mean(chisquareds))
    print("StdDev: ", scipy.std(chisquareds))
    print("Samples: ", len(chisquareds))
    print("Max: ", max(chisquareds))
    print("Min: ", min(chisquareds))
    bins = 20
    counts = [0 for i in range(bins)]
    for c in chisquareds:
        counts[int(c*bins)] += 1
    plot(counts)
#    plot(alignswith3)
