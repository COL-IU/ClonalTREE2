from predict import *
from math import *
import sys
from myutils import *
from graphviz import Digraph
from scipy import stats


def FIT_t(nus, time):
    start = 0
    for i in range(0, len(nus)):
        if 0 < nus[i] < 1:
            start = i
            break
    end = 0
    for i in range(len(nus) - 1, -1, -1):
        if 0 < nus[i] < 1:
            end = i
            break
    sub_nus = nus[start: end + 1]
    if len(sub_nus) <= 2:
        return zero, 0, 0, 0
    for nu in sub_nus:
        if nu == zero or nu == one:
            return zero, 0, 0, 0
    L = Decimal(len(sub_nus) - 1)
    Ys = []
    for i in range(1, len(sub_nus)):
        Y = (sub_nus[i] - sub_nus[i - 1]) / Decimal(sqrt(2 * sub_nus[i - 1] * (1 - sub_nus[i - 1]) * time))
        Ys.append(Y)
    Y_mean = sum(Ys) / L
    if Y_mean == zero: # Validate this decision later
        return zero, 0, 0, 0
    temp = zero
    for Y in Ys:
        temp = temp + ((Y - Y_mean)**2)
    S2 = temp / (L - 1)
    t = Y_mean / Decimal(sqrt(S2 / L))
    p = 2 * (1 - stats.t.cdf(float(abs(t)), df=int(L-1)))
    return t, L - 1, p, 1


if len(sys.argv) <= 2:
    sys.exit("Usage: python3 ClonalTREE2.py <prefix> <generation_time> <optional:gff>\n\n"
             "<prefix>:\t[String] Path/filename prefix for the input files and the output files.\n"
             "<generation_time>:\t[Int] Number of generations between each time point.\n\n"
             "Input files:\n"
             "<prefix>.vaf:\tInput file containing the variant allele frequencies matrix (F).\n"
             "<prefix>.rd:\tInput file containing the read depth matrix (R).\n"
             "<prefix>.var:\tInput file containing the variant names / loci.\n"
             "optional:file.gff:\tGFF3 file containing the gene annotation of the reference genome.\n\n"
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
gff_file = ""
time = int(sys.argv[2])

if len(sys.argv) == 4:
    gff_file = sys.argv[3]
    annotate = True
else:
    annotate = False

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

parents, score, variants1, removed_variants, num_times, running_time, removed_time_points, F1, R1, order = predict(F, variants, 1, R)

f = open(prefix + ".tree", "w")
write_parents(variants1, parents, f)
f.close()

my_F = add_founder(F1)
my_R = add_founder(R1)
# steps, arrival_times = get_step_structure(my_F)
C, F2 = get_c_no_fail(my_F, parents, order)

first_row = ["0"] + variants1

for_analysis = [0] * len(first_row)
f = open(prefix + ".C", "w")
C = [first_row] + C
C_t = list(map(list, zip(*C)))
C_ts = []
C_dfs = []
C_ps = []
for i in range(0, len(C_t)):
    t, df, p, fa = FIT_t(C_t[i][1:], time)
    if fa:
        if p < 0.1:
            for_analysis[i] = 1
        C_ts.append("{0:.3f}".format(t))
        C_dfs.append(str(df))
        C_ps.append("{0:.3f}".format(p))
    else:
        C_ts.append(" ")
        C_dfs.append(" ")
        C_ps.append(" ")
    # C_ts.append(" ")
    # C_dfs.append(" ")
    # C_ps.append(" ")
    # print(C_t[i], num_positive)
# C_t.sort(key=lambda x: x[-1], reverse=True)
C_sorted = list(map(list, zip(*C_t)))
f.write("\t".join(C_sorted[0]) + "\n")
write_dm(C_sorted[1:], f)
f.write("\n")
f.write("\t".join(C_ts))
f.write("\n")
f.write("\t".join(C_dfs))
f.write("\n")
f.write("\t".join(C_ps))
f.close()

f = open(prefix + ".F", "w")
F2 = [first_row] + F2
F_t = list(map(list, zip(*F2)))
F_ts = []
F_dfs = []
F_ps = []
for i in range(0, len(F_t)):
    t, df, p, fa = FIT_t(F_t[i][1:], time)
    if fa:
        F_ts.append("{0:.3f}".format(t))
        F_dfs.append(str(df))
        F_ps.append("{0:.3f}".format(p))
    else:
        F_ts.append(" ")
        F_dfs.append(" ")
        F_ps.append(" ")
    # print(C_t[i], num_positive)
# C_t.sort(key=lambda x: x[-1], reverse=True)
F_sorted = list(map(list, zip(*F_t)))
f.write("\t".join(F_sorted[0]) + "\n")
write_dm(F_sorted[1:], f)
f.write("\n")
f.write("\t".join(F_ts))
f.write("\n")
f.write("\t".join(F_dfs))
f.write("\n")
f.write("\t".join(F_ps))
f.close()

# f = open(prefix + ".F", "w")
# f.write("\t".join(first_row) + "\n")
# write_dm(my_F, f)
# f.close()

f = open(prefix + ".R", "w")
f.write("\t".join(first_row) + "\n")
write_dm(my_R, f)
f.close()

f = open(prefix + ".info", "w")
f.write("Score: " + str(score) + "\n")
f.write("Running Time: " + str(running_time) + "\n")
f.write("Removed Variants: " + str(removed_variants) + "\n")
f.write("Removed Time Points: " + str(removed_time_points) + "\n")
f.close()

genes = {}
if annotate:
    regions = []
    f = open(gff_file)
    for line in f:
        if line[0] == "#":
            continue
        words = line.split(None)
        if words[2] != "CDS":
            continue
        annotations = words[8].split(";")
        for annotation in annotations:
            if annotation[:5] == "gene=":
                regions.append([int(words[3]), int(words[4]), annotation[5:]])
                break
    f.close()
    for key in loci.keys():
        if key == "0":
            continue
        gene_names = ""
        variant_names = loci[key]
        variant_names = variant_names.split(",")
        for variant_name in variant_names:
            words = variant_name.split("@")
            position = int(words[1])
            for region in regions:
                if region[0] <= position <= region[1]:
                    gene_names = gene_names + region[2] + ","
                    break
        genes[key] = gene_names[0:-1]

dot = Digraph()
# if "0" in C_sorted[0][0:20]:
if for_analysis[0]:
    dot.node("0", "0 / Founder / 0\n" + "t=" + C_ts[0] + ",p=" + C_ps[0], style="filled", fillcolor="plum3")
else:
    dot.node("0", "0 / Founder / 0")
# for i in variants:
#     dot.node(str(i))
for i in range(0, len(variants1)):
    variant = variants1[i]
    # if variant in C_sorted[0][0:20]:
    if for_analysis[i + 1]:
        dot.attr('node', style="filled", fillcolor="plum3")
    else:
        dot.attr('node', style="filled", fillcolor="white")
    if annotate:
        if for_analysis[i + 1]:
            dot.node(str(variant), str(variant) + " / " + loci[str(variant)] + "\n"
                     + "t=" + C_ts[i + 1] + ",p=" + C_ps[i + 1] + "\n" + genes[str(variant)])
        else:
            dot.node(str(variant), str(variant) + " / " + loci[str(variant)] + "\n"
                     + genes[str(variant)])
    else:
        if for_analysis[i + 1]:
            dot.node(str(variant), str(variant) + " / " + loci[str(variant)] + "\n"
                     + "t=" + C_ts[i + 1] + ",p=" + C_ps[i + 1])
        else:
            dot.node(str(variant), str(variant) + " / " + loci[str(variant)])
for i in range(0, len(variants1)):
    to_node = str(variants1[i])
    if parents[i + 1] == 0:
        from_node = "0"
    else:
        from_node = str(variants1[parents[i + 1] - 1])
    dot.edge(from_node, to_node)
f = open(prefix + ".dot", "w")
f.write(dot.source)
f.close()
