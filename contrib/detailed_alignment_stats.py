#!/usr/bin/python3

import re
import sys
import scipy


alignment_line_re = re.compile(r"\[(.*)\]->\[(.*)\] s: (.*)")
frag_re = re.compile(r"(\d+):(\d+(\.\d+)?)")

if __name__ == "__main__":
    devs = []
    for line in open(sys.argv[1]):
        mobj = alignment_line_re.search(line.strip())
        if mobj:
            frags_l = frag_re.findall(mobj.group(1))
            frags_r = frag_re.findall(mobj.group(2))
            if len(frags_l) < 1 or len(frags_r) < 1:
                print("Unable to parse one or more frag subseqs in", mobj.groups()[:2])

            #if len(frags_l) > 1 or len(frags_r) > 1 : continue # uncomment to get only 1:1 alignments
            l_fragsum = sum((float(fstring[1]) for fstring in frags_l))
            r_fragsum = sum((float(fstring[1]) for fstring in frags_r))
            if (l_fragsum > 4):
                devs.append((l_fragsum - r_fragsum) / l_fragsum)

    print("Mean: ", scipy.mean(devs))
    print("StdDev: ", scipy.std(devs))
    print("Samples: ", len(devs))
    
