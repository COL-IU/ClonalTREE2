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


def cluster(F, variants, threshold):
    F_t = list(map(list, zip(*F)))
    clusters = [(variants[0], F_t[0])]
    i = 1
    while i < len(F_t):
        new_cluster = True
        for j in range(0, len(clusters)):
            (header, column) = clusters[j]
            if diff_below_cutoff(F_t[i], column, threshold):
                clusters[j] = (header + "_" + variants[i], average_column(F_t[i], column))
                new_cluster = False
                break
        if new_cluster:
            clusters.append((variants[i], F_t[i]))
        i = i + 1
    cluster_variants = []
    cluster_F = []
    for (cluster_variant, column) in clusters:
        cluster_variants.append(cluster_variant)
        cluster_F.append(column)
    return cluster_variants, list(map(list, zip(*cluster_F)))


if len(sys.argv) <= 2:
    sys.exit("cluster.py <threshold> <prefix>")

threshold = float(sys.argv[1])
arg1 = sys.argv[2] + ".vaf"
arg2 = sys.argv[2] + ".var"

f = open(arg1)
F, _ = read_F(f)
f.close()

f = open(arg2)
lines = f.readlines()
f.close()

variants = []
for line in lines:
    variants.append(line.strip())

(cluster_variants, cluster_F) = cluster(F, variants, threshold)

print(len(variants))
print(len(cluster_variants))

f = open(arg2 + "1", "w")
for cluster_variant in cluster_variants:
    f.write(cluster_variant + "\n")
f.close()

f = open(arg1 + "1", "w")
for i in range(0, len(cluster_variants)):
    f.write(str(i + 1) + "\t")
f.write("\n")
write_dm(cluster_F, f)
f.close()

