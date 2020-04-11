import sys
from decimal import *


def read_F(infile):
    F = []
    lines = infile.readlines()
    var = lines[0].split(None)
    for i in range(1, len(lines)):
        line = lines[i]
        words = line.split(None)
        F.append(list(map(Decimal, words)))
    return F, var


def write_dm(mat, out):
    m = len(mat)
    n = len(mat[0])
    for i in range(0, m):
        for j in range(0, n):
            out.write("{0:.5f}".format(mat[i][j]) + "\t")
        out.write("\n")


if len(sys.argv) <= 1:
    sys.exit("get_non_syn.py <prefix>")

arg1 = sys.argv[1] + ".vaf"
arg2 = sys.argv[1] + ".var"
arg3 = sys.argv[1] + ".putations.K12MG1655_Chr1.detail"
arg4 = sys.argv[1] + ".rd"
out_prefix = sys.argv[1]

f = open(arg1)
F, _ = read_F(f)
f.close()

f = open(arg2)
lines = f.readlines()
f.close()

f = open(arg4)
R, _ = read_F(f)
f.close()

variants = []
for line in lines:
    variants.append(line.strip())

anno = {}
f = open(arg3)
lines = f.readlines()
for i in range(1, len(lines)):
    line = lines[i]
    words = line.split(None)
    anno[words[2] + "@" + words[3]] = words[16]
f.close()

F_t = list(map(list, zip(*F)))
R_t = list(map(list, zip(*R)))
variants_new = []
F_new_t = []
R_new_t = []

for i in range(0, len(variants)):
    variant = variants[i]
    if variant in anno.keys():
        if anno[variant] == "N-Syn":
            variants_new.append(variant)
            F_new_t.append(F_t[i])
            R_new_t.append(R_t[i])

F_new = list(map(list, zip(*F_new_t)))
R_new = list(map(list, zip(*R_new_t)))

f = open(out_prefix + ".ns.var", "w")
for variant in variants_new:
    f.write(variant + "\n")
f.close()

f = open(out_prefix + ".ns.vaf", "w")
for i in range(0, len(variants_new)):
    f.write(str(i + 1) + "\t")
f.write("\n")
write_dm(F_new, f)
f.close()

f = open(out_prefix + ".ns.rd", "w")
for i in range(0, len(variants_new)):
    f.write(str(i + 1) + "\t")
f.write("\n")
write_dm(R_new, f)
f.close()
