from predict import *
from math import *
import sys
from myutils import *
from graphviz import Digraph

if len(sys.argv) <= 1:
    sys.exit("Usage: python3 ClonalTREE.py <prefix>\n\n"
             "<prefix>:\t[String] Path/filename prefix for the input files and the output files.\n\n"
             "Input files:\n"
             "<prefix>.vaf:\tInput file containing the variant allele frequencies matrix (F).\n"
             "<prefix>.rd:\tInput file containing the read depth matrix (R).\n"
             "<prefix>.var:\tInput file containing the variant names / loci.\n\n"
             "Output files:\n"
             "<prefix>.tree:\tList of each node and their corresponding ancestor.\n"
             "<prefix>.dot:\tTree in dot format to visualize using GraphViz.\n"
             "<prefix>.C:\tTable containing the clonal frequency matrix (C).\n"
             "<prefix>.info:\tA few added information regarding the prediction.\n"
             )

vaf_file = sys.argv[1] + ".vaf"
rd_file = sys.argv[1] + ".rd"
var_file = sys.argv[1] + ".var"
prefix = sys.argv[1]

# arg1 = "pop125.vaf"
# arg2 = "pop125.rd"
# arg3 = "pop125.var"
# arg4 = "pop125"

f = open(vaf_file)
F, variants = read_F(f)
f.close()

f = open(rd_file)
R, _ = read_F(f)
f.close()

loci = {"0": "Founder"}
f = open(var_file)
lines = f.readlines()
for i in range(0, len(lines)):
    line = lines[i]
    loci[str(i + 1)] = line.strip()
f.close()

parents, score, variants, removed_variants, num_times, running_time, removed_time_points, F = predict(F, variants, 1, R)

f = open(prefix + ".tree", "w")
write_parents(variants, parents, f)
f.close()

my_F = add_founder(F)
C = get_c_no_fail(my_F, parents)

f = open(prefix + ".C", "w")
C_t = list(map(list, zip(*C)))
C_t[0] = ["0"] + C_t[0]
for i in range(1, len(C_t)):
    C_t[i] = [variants[i-1]] + C_t[i]
C_t.sort(key=lambda x: x[-1], reverse=True)
C = list(map(list, zip(*C_t)))
f.write("\t".join(C[0]) + "\n")
write_dm(C[1:], f)
f.close()

f = open(prefix + ".info", "w")
f.write("Score: " + str(score) + "\n")
f.write("Running Time: " + str(running_time) + "\n")
f.write("Removed Variants: " + str(removed_variants) + "\n")
f.write("Removed Time Points: " + str(removed_time_points) + "\n")
f.close()

dot = Digraph()
dot.node("0")
# for i in variants:
#     dot.node(str(i))
for i in variants:
    dot.node(str(i), loci[str(i)])
for i in range(0, len(variants)):
    to_node = str(variants[i])
    if parents[i + 1] == 0:
        from_node = "0"
    else:
        from_node = str(variants[parents[i + 1] - 1])
    dot.edge(from_node, to_node)
f = open(prefix + ".dot", "w")
f.write(dot.source)
f.close()
