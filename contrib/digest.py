#!/usr/bin/python3
import sys
import Bio
import Bio.SeqIO
import Bio.Restriction

if len(sys.argv) < 4:
    print("Usage:", sys.argv[0], "<reference.fa> <enzyme> <silico_map.valuev>")
    sys.exit(1)


for record in Bio.SeqIO.parse(open(sys.argv[1]), "fasta"):
    ref_seq = record.seq

enzyme_name = sys.argv[2]

lengths = [len(s) for s in getattr(Bio.Restriction, enzyme_name).catalyze(ref_seq)]
silico = open(sys.argv[3], "w")
silico.write(sys.argv[1] + "\n")
silico.write("\t" + enzyme_name + "\t" + enzyme_name + "\t" + "\t".join([str(float(x)/1000.0) for x in lengths]))
silico.write("\n")
silico.close()
