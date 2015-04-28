import re
import sys
import scipy

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

if __name__ == "__main__":
    matches = []
    num_frags = []
    nontriv_matches = []
    for line in open(sys.argv[1]):
        if align_pair_re.search(line):
            matches.append(0)
            nontriv_matches.append(0)
        mobj = alignment_line_re.search(line)
        if mobj:
            frags_l = frag_re.findall(mobj.group(1))
            frags_r = frag_re.findall(mobj.group(2))
            incr(num_frags, len(frags_l))
            incr(num_frags,len(frags_r))
            if len(frags_l) > 1:
                nontriv_matches[-1] += 1

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

