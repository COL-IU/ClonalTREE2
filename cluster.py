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


def get_max_diff(l1, l2):
    return max(map(lambda x, y: abs(x - y), l1, l2))


def diff_below_cutoff(l1, l2, cutoff):
    if len(l1) != len(l2):
        return False
    for i in range(0, len(l1)):
        if abs(l1[i] - l2[i]) > cutoff:
            return False
    return True


def average_column(l1, l2):
    # Assume equal length
    return list(map(lambda x, y: (x + y) / 2, l1, l2))


def cluster(F, variants, threshold, R):
    F_t = list(map(list, zip(*F)))
    R_t = list(map(list, zip(*R)))
    clusters = [(variants[0], F_t[0], R_t[0])]
    i = 1
    while i < len(F_t):
        new_cluster = True
        for j in range(0, len(clusters)):
            (header, column, r_column) = clusters[j]
            if diff_below_cutoff(F_t[i], column, threshold):
                clusters[j] = (header + "," + variants[i], average_column(F_t[i], column), average_column(R_t[i], r_column))
                new_cluster = False
                break
        if new_cluster:
            clusters.append((variants[i], F_t[i], R_t[i]))
        i = i + 1
    cluster_variants = []
    cluster_F = []
    cluster_R = []
    for (cluster_variant, column, r_column) in clusters:
        cluster_variants.append(cluster_variant)
        cluster_F.append(column)
        cluster_R.append(r_column)
    return cluster_variants, list(map(list, zip(*cluster_F))), list(map(list, zip(*cluster_R)))


if len(sys.argv) <= 2:
    sys.exit("cluster.py <threshold> <prefix>")

threshold = float(sys.argv[1])
arg1 = sys.argv[2] + ".vaf"
arg2 = sys.argv[2] + ".var"
arg3 = sys.argv[2] + ".rd"

f = open(arg1)
F, _ = read_F(f)
f.close()

f = open(arg2)
lines = f.readlines()
f.close()

f = open(arg3)
R, _ = read_F(f)
f.close()

variants = []
for line in lines:
    variants.append(line.strip())

(cluster_variants, cluster_F, cluster_R) = cluster(F, variants, threshold, R)

# print(len(variants))
# print(len(cluster_variants))

out1 = sys.argv[2] + ".cs.vaf"
out2 = sys.argv[2] + ".cs.var"
out3 = sys.argv[2] + ".cs.rd"

f = open(out2, "w")
for cluster_variant in cluster_variants:
    f.write(cluster_variant + "\n")
f.close()

f = open(out1, "w")
for i in range(0, len(cluster_variants)):
    f.write(str(i + 1) + "\t")
f.write("\n")
write_dm(cluster_F, f)
f.close()

f = open(out3, "w")
for i in range(0, len(cluster_variants)):
    f.write(str(i + 1) + "\t")
f.write("\n")
write_dm(cluster_R, f)
f.close()
