from predict import *
from math import *
import sys
from myutils import *
from graphviz import Digraph

if len(sys.argv) <= 3:
    sys.exit("Usage: python3 ClonalTREE.py <VAF file> <rd file> <out prefix>\n\n"
             "VAF file:\t[String] Input file containing the variant allele frequencies matrix (F).\n"
             "RD file:\t[String] Input file containing the read depth matrix (R).\n"
             "out prefix:\t[String] File path to prefix all output files.\n")

arg1 = sys.argv[1]
arg2 = sys.argv[2]
arg3 = sys.argv[3]

# arg1 = "pop125.vaf"
# arg2 = "pop125.rd"
# arg3 = "pop125"

f1 = open(arg1)
F, variants = read_F(f1)
f1.close()

f2 = open(arg2)
R, _ = read_F(f2)
f2.close()

out_prefix = arg3

parents, score, variants, removed_variants, num_times, running_time, removed_time_points = predict(F, variants, 1, R)

o1 = open(out_prefix + ".tree", "w")
write_parents(variants, parents, o1)
o1.close()

o2 = open(out_prefix + ".info", "w")
o2.write("Score: " + str(score) + "\n")
o2.write("Running Time: " + str(running_time) + "\n")
o2.write("Removed Variants: " + str(removed_variants) + "\n")
o2.write("Removed Time Points: " + str(removed_time_points) + "\n")
o2.close()

dot = Digraph()
for i in range(0, len(variants)):
    dot.node(str(i))
for i in range(0, len(variants)):
    to_node = str(variants[i])
    if parents[i + 1] == 0:
        from_node = "0"
    else:
        from_node = str(variants[parents[i + 1] - 1])
    dot.edge(from_node, to_node)
o3 = open(out_prefix + ".dot", "w")
o3.write(dot.source)
o3.close()
