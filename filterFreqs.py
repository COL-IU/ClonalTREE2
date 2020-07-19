import sys
from decimal import *
from copy import *


zero = Decimal('0.0')


def filter_founder_variants(F, R, variants):
    F_temp = deepcopy(F)
    R_temp = deepcopy(R)
    F_temp = list(map(list, zip(*F_temp)))
    R_temp = list(map(list, zip(*R_temp)))
    F_out = []
    R_out = []
    founder_F_out = []
    founder_R_out = []
    variants_out = []
    founder_var = []
    # print(len(variants))
    for i in range(0, len(F_temp)):
        num_ones1 = 0
        for j in range(0, len(F_temp[i])):
            if F_temp[i][j] > 0.95:
                num_ones1 += 1
        if num_ones1 < 4:
            F_out.append(F_temp[i])
            R_out.append(R_temp[i])
            variants_out.append(variants[i])
        else:
            founder_F_out.append(F_temp[i])
            founder_R_out.append(R_temp[i])
            founder_var.append(variants[i])
    # print(len(variants_out))
    # print(len(founder_var))
    F_out = list(map(list, zip(*F_out)))
    R_out = list(map(list, zip(*R_out)))
    founder_F_out = list(map(list, zip(*founder_F_out)))
    founder_R_out = list(map(list, zip(*founder_R_out)))
    return F_out, R_out, variants_out, founder_var, founder_F_out, founder_R_out


def reorder(l, order):
    return [l[i] for i in order]


def rearrange(F, R, variants, remove_redundancy=False):
    F_out = deepcopy(F)
    R_temp = deepcopy(R)
    R_out = []
    m = len(F_out)
    if m == 0:
        return F, R, variants
    if len(F[0]) == 0:
        return F, R, variants
    # Make F_out a diagonal (step) matrix
    F_out = list(map(list, zip(*F_out)))
    R_temp = list(map(list, zip(*R_temp)))
    order = []
    i = 0
    for row in F_out:
        for i in range(0, len(row)):
            if row[i] > zero:
                break
        order.append(i)
    indices = list(range(0, len(F_out)))
    [order, indices, F_out] = list(zip(*sorted(zip(order, indices, F_out), key=lambda x: x[0])))
    for index in indices:
        R_out.append(R_temp[index])
    variants_out = reorder(variants, indices)

    non_founder_F = []
    non_founder_var = []
    non_founder_R = []
    founder_var = []

    for i in range(0, len(variants_out)):
        row = F_out[i]
        is_founder = True
        for j in range(0, len(row)):
            if row[j] < Decimal('0.9'):
                is_founder = False
                break
        if is_founder:
            founder_var.append(variants_out[i])
        else:
            non_founder_F.append(row)
            non_founder_var.append(variants_out[i])
            non_founder_R.append(R_out[i])

    non_founder_F = list(map(list, zip(*non_founder_F)))
    non_founder_R = list(map(list, zip(*non_founder_R)))
    # non_founder_F = list(map(list, zip(*F_out)))
    # non_founder_R = list(map(list, zip(*R_out)))

    if remove_redundancy:
        # Remove unnecessary time points
        z = 0
        z_prev = -1
        new_F = []
        new_R = []
        for j in range(0, len(non_founder_F)):
            row = non_founder_F[j]
            for z in range(len(row)-1, -2, -1):
                if row[z] > zero:
                    break
            if z - z_prev > 0:
                new_F.append(non_founder_F[j])
                new_R.append(non_founder_R[j])
            z_prev = z
        non_founder_F = new_F
        non_founder_R = new_R

    return non_founder_F, non_founder_R, non_founder_var, founder_var


def is_valid(t):
    if len(t) < 2:
        return False
    ts = sorted(t)
    # if max(ts) != 6:
    #     return False
    for i in range(1, len(ts)):
        if ts[i] - ts[i-1] > 1:
            return False
    return True


def write_dm(mat, out):
    m = len(mat)
    n = len(mat[0])
    for i in range(0, m):
        for j in range(0, n):
            out.write("{0:.5f}".format(mat[i][j]) + "\t")
        out.write("\n")


if len(sys.argv) < 3:
    sys.exit("Usage: filterFreqs.py <freqs> <outprefix>")

f = open(sys.argv[1])
lines = f.readlines()
f.close()

times = {}
subs = {}
invalids = []
max_time = 0

for line in lines:
    words = line.split(None)
    locus = words[1]
    time = int(words[0][1:])
    if time > max_time:
        max_time = time
    if locus in times:
        times[locus].append(time)
    else:
        times[locus] = [time]

    if locus in subs.keys():
        if subs[locus] != (words[3], words[4]):
            invalids.append(locus)
            # print(line)
    else:
        subs[locus] = (words[3], words[4])

filtered = []

for locus in sorted(times.keys()):
    if is_valid(times[locus]):
        filtered.append(locus)

f1 = open(sys.argv[2]+".filtered.freqs", "w")
F = [[zero for _ in range(len(filtered))] for _ in range(max_time)]
R = [[zero for _ in range(len(filtered))] for _ in range(max_time)]
for line in lines:
    words = line.split(None)
    locus = words[1]
    time = int(words[0][1:])-1
    if locus in filtered:
        col = filtered.index(locus)
        F[time][col] = Decimal(words[6])
        R[time][col] = Decimal(words[2])
        f1.write(line.strip() + "\n")
f1.close()

F, R, filtered, founder_variants, founder_F, founder_R = filter_founder_variants(F, R, filtered)
F_new, R_new, var_new, more_founder = rearrange(F, R, filtered, False)

founder_variants = founder_variants + more_founder

f2 = open(sys.argv[2]+".vaf", "w")
if len(F_new) > 0:
    n = len(F_new[0])
    f2.write("\t".join(list(map(str, range(1,n+1)))))
    f2.write("\n")
    write_dm(F_new, f2)
f2.close()

f3 = open(sys.argv[2]+".var", "w")
for var in var_new:
    f3.write(var + "\n")
f3.close()

f4 = open(sys.argv[2]+".rd", "w")
if len(R_new) > 0:
    n = len(R_new[0])
    f4.write("\t".join(list(map(str, range(1,n+1)))))
    f4.write("\n")
    write_dm(R_new, f4)
f4.close()

f5 = open(sys.argv[2]+".putations", "w")
variants = []
for variant in var_new:
    if variant not in invalids:
        temp = variant.split("@")
        contig = temp[0]
        locus = int(temp[1])
        variants.append([contig, locus])
for [contig, locus] in sorted(variants, key=lambda x: (x[0], x[1])):
    key = contig + "@" + str(locus)
    sub = subs[key]
    f5.write(contig + " " + str(locus) + " " + sub[0] + " " + sub[1] + " " + sub[0] + "\n")
f5.close()

f6 = open(sys.argv[2]+".fv.putations", "w")
variants = []
for variant in founder_variants:
    if variant not in invalids:
        temp = variant.split("@")
        contig = temp[0]
        locus = int(temp[1])
        variants.append([contig, locus])
for [contig, locus] in sorted(variants, key=lambda x: (x[0], x[1])):
    key = contig + "@" + str(locus)
    sub = subs[key]
    f6.write(contig + " " + str(locus) + " " + sub[0] + " " + sub[1] + " " + sub[0] + "\n")
f6.close()

f = open(sys.argv[2]+".fv.vaf", "w")
if len(founder_F) > 0:
    write_dm(founder_F, f)
f.close()

f = open(sys.argv[2]+".fv.R", "w")
if len(founder_R) > 0:
    write_dm(founder_R, f)
f.close()
