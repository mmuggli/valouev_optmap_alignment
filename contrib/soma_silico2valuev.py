#!/usr/bin/python

# makes a valouev style optical map from soma's in-silico digest


import sys


if len(sys.argv) < 2:
    print "Usage: " + sys.argv[0] + " <soma format in silico digested contigs file>      >     valouev_frag_seqs_file"
    sys.exit(1)

def proc(line):
    parts = line.split(' ')
    newparts = [ str(int(part) / 1000.0 ) for part in parts] + [length / 1000.0]
    newdiffs = [str(j - i) for i,j in [(float(x), float(y)) for x,y in zip(newparts, newparts[1:])]]
    return " ".join(newdiffs)

for ln, line in enumerate(open(sys.argv[1])):
    if ln % 2 == 1: 
        print "enzyme enzyme",proc(line.strip())
    else:
        if ln % 2 == 0: length = int(line.strip().split()[1])
        print line, 
    if ln % 2 == 1: print
