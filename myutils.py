from sys import *
from decimal import *
from itertools import *
from copy import *
from math import *
from time import *
from datetime import *
import numpy as np
from datetime import datetime
import networkx as nx

one = Decimal('1.0')
zero = Decimal('0.0')
low = Decimal('0.00001')
high = Decimal('Inf')
FAIL_THRESHOLD = Decimal('0.0')
high_read_depth = maxsize


def rd_to_sd(F, R):
    m = len(R)
    if m == 0:
        S = []
    else:
        n = len(R[0])
        S = [[Decimal('0.0') for _ in range(n)] for _ in range(m)]
        for i in range(0, m):
            for j in range(0, n):
                if R[i][j] == zero:
                    S[i][j] = zero
                else:
                    if F[i][j] == one:
                        S[i][j] = Decimal(sqrt((low * (1 - low)) / R[i][j]))
                    else:
                        S[i][j] = Decimal(sqrt((F[i][j] * (1 - F[i][j])) / R[i][j]))
    return S


def compare_columns_undirected(c1, c2, ind1, ind2):
    diff1 = 0
    count1 = 0
    diff2 = 0
    count2 = 0
    for i in range(0, len(c1)):
        if c1[i] >= c2[i]:
            diff1 = diff1 + c1[i] - c2[i]
            count1 = count1 + 1
        if c2[i] >= c1[i]:
            diff2 = diff2 + c2[i] - c1[i]
            count2 = count2 + 1
    if count1 == len(c1):
        return [(ind1, ind2, 0)]
    if count2 == len(c1):
        return [(ind2, ind1, 0)]
    if diff1 > diff2:
        return [(ind1, ind2, diff2)]
    else:
        return [(ind2, ind1, diff1)]


def compare_columns_directed(c1, c2, ind1, ind2):
    diff1 = 0
    count1 = 0
    diff2 = 0
    count2 = 0
    for i in range(0, len(c1)):
        if c1[i] >= c2[i]:
            diff1 = diff1 + c1[i] - c2[i]
            count1 = count1 + 1
        if c2[i] >= c1[i]:
            diff2 = diff2 + c2[i] - c1[i]
            count2 = count2 + 1
    if count1 == len(c1):
        return [(ind1, ind2, 0)]
    if count2 == len(c1):
        return [(ind2, ind1, 0)]
    return [(ind1, ind2, diff2), (ind2, ind1, diff1)]


def get_partial_order(F, column_range):
    # print_dm(F)
    F_t = list(map(list, zip(*F)))
    G = {}
    vertices = set()
    weighted_edges_undirected = []
    weighted_edges_directed = []
    edges = []
    for i in range(column_range[0], column_range[1]):
        for j in range(i + 1, column_range[1] + 1):
            # weighted_edge_undirected = compare_columns_undirected(F_t[i], F_t[j], i, j)
            # weighted_edges_undirected = weighted_edges_undirected + weighted_edge_undirected
            # edges.append((weighted_edge_undirected[0][0], weighted_edge_undirected[0][1]))
            weighted_edge_directed = compare_columns_directed(F_t[i], F_t[j], i, j)
            weighted_edges_directed = weighted_edges_directed + weighted_edge_directed

    # UG = nx.Graph()
    # UG.add_weighted_edges_from(weighted_edges_undirected)
    # mst = nx.minimum_spanning_tree(UG)
    # print(weighted_edges_undirected)
    # print(nx.get_edge_attributes(mst, "weight"))

    DG = nx.DiGraph()
    DG.add_weighted_edges_from(weighted_edges_directed)
    msa = nx.minimum_spanning_arborescence(DG, attr='weight', default=1)
    # print(weighted_edges_directed)
    # print(nx.get_edge_attributes(msa, "weight"))

    # mstEdges = list(map(sorted, mst.edges()))
    # for edge in edges:
    #     if sorted(edge) in mstEdges:
    #         start = edge[0]
    #         end = edge[1]
    #         if start in G.keys():
    #             G[start].add(end)
    #         else:
    #             G[start] = {end}
    #         vertices.add(start)
    #         vertices.add(end)
    # for vertex in vertices:
    #     if vertex not in G.keys():
    #         G[vertex] = set()

    for edge in msa.edges():
        start = edge[0]
        end = edge[1]
        if start in G.keys():
            G[start].add(end)
        else:
            G[start] = {end}
        vertices.add(start)
        vertices.add(end)
    for vertex in vertices:
        if vertex not in G.keys():
            G[vertex] = set()

    return G, list(vertices)


def topological_sorts_util(G, vertices, visited, in_degree, stack, out):
    flag = False
    for i in vertices:
        if not visited[i] and in_degree[i] == 0:
            visited[i] = True
            stack.append(i)
            for adjacent in G[i]:
                in_degree[adjacent] = in_degree[adjacent] - 1
            out = topological_sorts_util(G, vertices, visited, in_degree, stack, out)
            visited[i] = False
            stack.pop()
            for adjacent in G[i]:
                in_degree[adjacent] = in_degree[adjacent] + 1
            flag = True
    if not flag:
        out.append(stack[:])
        l = len(out)
        if l % 1000000 == 0:
            print(format(l, "10.2E"), datetime.now(), flush=True)
    return out


def topological_sorts(G, vertices):
    out = []
    visited = {}
    for i in vertices:
        visited[i] = False
    in_degrees = {}
    for i in vertices:
        in_degrees[i] = 0
    for i in vertices:
        for var in G[i]:
            in_degrees[var] = in_degrees[var] + 1
    stack = []
    out = topological_sorts_util(G, vertices, visited, in_degrees, stack, out)
    return out


def topological_valid_permutations(F, step_unit):
    step_start = step_unit[0]
    step_end = step_unit[-1]
    if step_start == step_end:
        return [[step_start]], F
    G, vertices = get_partial_order(F, (step_start, step_end))
    return topological_sorts(G, vertices)


def topological_sort_util(G, v, visited, stack):
    visited[v] = True
    for i in G[v]:
        if not visited[i]:
            topological_sort_util(G, i, visited, stack)
    stack.insert(0, v)


def topological_sort(G, vertices):
    visited = {}
    for i in vertices:
        visited[i] = False
    stack = []
    for i in vertices:
        if not visited[i]:
            topological_sort_util(G, i, visited, stack)
    return stack


def remove_first_rowcol(mat):
    out = []
    for i in range(1, len(mat)):
        out.append(mat[i][1:])
    return out


def add_founder_reads(R):
    my_R = [[zero for _ in range(len(R[0]) + 1)]]
    my_R[0][0] = high_read_depth
    for row in R:
        my_R.append([high_read_depth] + row)
    return my_R


def add_founder_penalty(S):
    my_S = [[zero for _ in range(len(S[0]) + 1)]]
    my_S[0][0] = zero
    for row in S:
        my_S.append([zero] + row)
    return my_S


def add_founder(F):
    my_F = [[zero for _ in range(len(F[0]) + 1)]]
    my_F[0][0] = one
    for row in F:
        my_F.append([one] + row)
    return my_F


def ancestry_graph(F):
    m = len(F)
    n = len(F[0])
    G = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(0, n):
        edges = set(list(range(n)))
        for j in range(0, m):
            temp = set()
            for k in range(0, n):
                if F[j][k] <= F[j][i]:
                    temp.add(k)
            edges = edges & temp
        for e in edges:
            G[e][i] = 1
    return G


def enum_spanning_trees(G_i):
    G = deepcopy(G_i)
    n = len(G)
    for i in range(0, n):
        G[i][i] = 0
    ones = []
    for i in range(1, n):
        indices = [j for j, x in enumerate(G[i]) if x == 1]
        ones.append(indices)
    num_trees = 1
    for item in ones:
        num_trees = num_trees * len(item)
    stdout.write(" " + str(num_trees) + " trees\n")
    stdout.flush()
    sTrees = []
    i = 0
    for element in product(*ones):
        # if i % 1000000 == 0:
        #     stdout.write(" i1 = " + str(i/1000000) + " M; " + str(datetime.now()) + "\n")
        #     stdout.flush()
        # i += 1
        sTrees.append(list_to_tree([0] + list(element)))
    return sTrees


def calc_dist(l):
    dist = []
    s = Decimal(sum(l))
    for i in l:
        dist.append(Decimal(i) / s)
    return dist


def list_to_tree(parents):
    T = [[0] * len(parents)]
    T[0][0] = 1
    for i in range(1, len(parents)):
        parent = parents[i]
        T.append(T[parent][:])
        T[i][i] = 1
    return T


def list_to_clones(parents, variants):
    clones = [[]]
    for i in range(1, len(parents)):
        parent = parents[i]
        clones.append(clones[parent][:])
        clones[i].append(variants[i - 1])
    return clones[1:]


def tree_to_list(T):
    n = len(T)
    parents = [0] * n
    children = [[] for _ in range(n)]
    for i in range(1, n):
        row = T[i][:]
        row[i] = 0
        for j in range(0, n):
            if T[j] == row:
                parents[i] = j
                children[j].append(i)
    return parents, children


def get_c_penalty(F, S, T, fail_threshold):
    # Assumes square F
    n = len(T)
    parents, children = tree_to_list(T)
    C = [[zero for _ in range(n)] for _ in range(n)]
    C[0][0] = one
    fail = False
    for i in range(1, n):
        for j in range(0, n):
            c = children[j]
            cs = 0
            cs_sds = 0
            for child in c:
                cs += F[i][child]
                cs_sds += (S[i][child] ** 2)
            C[i][j] = F[i][j] - cs
            if S[i][j] > zero:
                C[i][j] = (F[i][j] - cs) + (np.exp(abs(F[i][j] - cs) / Decimal(sqrt((S[i][j] ** 2) + cs_sds))))
            if C[i][j] < fail_threshold:
                fail = True
    if fail:
        p = -1
    else:
        p = one
        for i in range(1, n):
            p = p * C[i - 1][parents[i]]
    return C, p


def get_c(F, T, fail_threshold):
    # Assumes square F
    n = len(T)
    parents, children = tree_to_list(T)
    C = [[zero for _ in range(n)] for _ in range(n)]
    C[0][0] = one
    fail = False
    for i in range(1, n):
        for j in range(0, n):
            c = children[j]
            cs = 0
            for child in c:
                cs += F[i][child]
            C[i][j] = F[i][j] - cs
            if C[i][j] < fail_threshold:
                fail = True
    if fail:
        p = -1
    else:
        p = one
        for i in range(1, n):
            p = p * C[i - 1][parents[i]]
    return C, p


def get_c_no_fail(F, parents, order):
    if not F:
        return [], []
    F1 = deepcopy(F)
    m = len(F1)
    n = len(F1[0])
    children_dict = parents_to_children(parents, len(parents))
    C = [[zero for _ in range(n)] for _ in range(m)]

    # Assume that F contains the founder column and row
    C[0][0] = one
    for i in range(1, m):
        for j in order:
            children = children_dict[j]
            cs = 0
            for child in children:
                cs += F1[i][child]
            C[i][j] = F1[i][j] - cs
            if C[i][j] < zero:
                C[i][j] = zero
                # F1[i][j] = (F1[i][j] + cs) / 2
                # if F1[i][j] > one:
                #     F1[i][j] = one
                # stdout.write("{0:.5f}".format(F1[i][j]) + "\t")
                for child in children:
                    # if (F1[i][child] / cs) * F1[i][j] == zero and F1[i][child] != zero:
                    #     print(i, j, child, F1[i][j], F1[i][child])
                    F1[i][child] = (F1[i][child] / cs) * F1[i][j]
    # for i in range(1, m):
    #     for j in order:
    #         children = children_dict[j]
    #         cs = 0
    #         for child in children:
    #             cs += F1[i][child]
    #         C[i][j] = F1[i][j] - cs

    # for i in range(1, m):
    #     s = sum(C[i])
    #     for j in range(0, n):
    #         C[i][j] = C[i][j] / s
    return C, F1


def get_p(F, parents, fail_threshold):
    tree = list_to_tree(parents)
    _, p = get_c(F, tree, fail_threshold)
    return p


def get_f(C, T):
    m = len(C)
    if m == 0:
        return []
    n = len(C[0])
    F = [[Decimal('0.0') for _ in range(n)] for _ in range(m)]
    for i in range(0, m):
        for j in range(0, n):
            temp = Decimal('0.0')
            for k in range(0, n):
                temp += (C[i][k] * T[k][j])
            F[i][j] = temp
    return F


def get_f_from_parents(C, parents):
    T = list_to_tree(parents)
    F = get_f(C, T)
    return F


def write_im(mat, out):
    m = len(mat)
    n = len(mat[0])
    for i in range(0, m):
        for j in range(0, n):
            out.write(str(mat[i][j]) + "\t")
        out.write("\n")


def write_dm(mat, out):
    m = len(mat)
    n = len(mat[0])
    for i in range(0, m):
        for j in range(0, n):
            out.write("{0:.5f}".format(mat[i][j]) + "\t")
        out.write("\n")


def read_F(infile):
    F = []
    lines = infile.readlines()
    var = lines[0].split(None)
    for i in range(1, len(lines)):
        line = lines[i]
        words = line.split(None)
        F.append(list(map(Decimal, words)))
    return F, var


def read_dm(infile):
    out = []
    lines = infile.readlines()
    for line in lines:
        words = line.split(None)
        out.append(list(map(Decimal, words)))
    return out


def print_dm(mat):
    # Assumes square mat
    m = len(mat)
    n = len(mat[0])
    for i in range(0, m):
        for j in range(0, n):
            stdout.write("{0:.5f}".format(mat[i][j]) + "\t")
        stdout.write("\n")


def print_dl(l):
    for item in l:
        stdout.write("{0:.3f}".format(item) + "\t")
    stdout.write("\n")


def dl_to_str(l):
    out = []
    for item in l:
        out.append("{0:.3f}".format(item))
    return "\t".join(out)


def print_tree(T):
    n = len(T)
    for i in range(0, n):
        for j in range(0, n):
            stdout.write(str(T[i][j]) + " ")
        stdout.write("\n")


def get_children(parents, parent):
    # print(parents, parent)
    children = []
    for i in range(1, len(parents)):
        if parents[i] == parent:
            children.append(i)
    # print(children)
    return children


def parents_to_children(parents, num_clones):
    children = {}
    for clone in range(0, num_clones):
        children[clone] = []
    for clone in range(1, num_clones):
        if parents[clone] != -1:
            children[parents[clone]].append(clone)
    return children


def children_to_parents(children, num_clones):
    parents = {}
    for clone in range(1, num_clones):
        parents[clone] = -1
    for clone in range(0, num_clones):
        cur_children = children[clone]
        for child in cur_children:
            parents[child] = clone
    return parents


def read_variants_set(f):
    out = set()
    for line in f:
        out.add(int(line.strip()))
    return out


def read_variants_list(f):
    out = []
    for line in f:
        out.append(int(line.strip()))
    return out


def reorder(l, order):
    return [l[i] for i in order]


def rearrange_rows(mat1, variants1, variants2):
    mat2 = []
    for var in variants2:
        mat2.append(mat1[variants1.index(var)])
    return mat2


def remove_redundant_time_points(F, R):
    selected_times = []
    removed_times = []
    z = 0
    z_prev = -1
    new_F = []
    new_R = []
    for j in range(0, len(F)):
        row = F[j]
        for z in range(len(row) - 1, -2, -1):
            if row[z] > zero:
                break
        if z - z_prev > 0 and not all([v == zero for v in row]):
            new_F.append(F[j])
            new_R.append(R[j])
            selected_times.append(j)
        else:
            removed_times.append(j)
        z_prev = z
    return new_F, new_R, selected_times, removed_times


def rearrange_penalty(F, R, variants):
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

    F_out = list(map(list, zip(*F_out)))
    R_out = list(map(list, zip(*R_out)))

    return F_out, R_out, variants_out, list(indices)


def rearrange(F, variants):
    F_out = deepcopy(F)
    m = len(F_out)
    if m == 0:
        return F, variants
    if len(F[0]) == 0:
        return F, variants
    # Make F_out a diagonal (step) matrix
    F_out = list(map(list, zip(*F_out)))
    order = []
    i = 0
    for row in F_out:
        for i in range(0, len(row)):
            if row[i] > zero:
                break
        order.append(i)
    indices = list(range(0, len(F_out)))
    [order, indices, F_out] = list(zip(*sorted(zip(order, indices, F_out), key=lambda x: x[0])))
    variants_out = reorder(variants, indices)

    F_out = list(map(list, zip(*F_out)))

    return F_out, variants_out, list(indices)


def get_step_structure(F):
    # Assumes step structured F
    # Also assumes that redundant time points are removed
    m = len(F)
    if m == 0:
        return []
    n = len(F[0])
    multi_spawns = []
    z_prev = -1
    for i in range(m):
        z = n - 1
        for z in range(n - 1, -1, -1):
            if F[i][z] > 0:
                break
        multi_spawns.append(range(z_prev + 1, z + 1))
        z_prev = z
    steps = list(map(list, multi_spawns))
    arrival_times = []
    i = 0
    for step in steps:
        arrival_times = arrival_times + ([i] * len(step))
        i = i + 1
    return steps, arrival_times


def correct_steps_for_known_founders(steps, k):
    out_steps = [[0]]
    step = steps[1]
    if k == 0 or k == step[-1]:
        return steps
    if step[-1] < k:
        exit("Invalid input for known founders. The VAF for founder variants in first time point should be non-zero.\n")
    out_steps.append(list(range(1, k + 1)))
    out_steps.append(list(range(k + 1, step[-1] + 1)))
    for i in range(2, len(steps)):
        out_steps.append(steps[i])
    return out_steps


def squarify(F, step_structure):
    # Assumes F's step structure matches with step_structure
    if not F:
        return []
    m = len(F)
    n = len(F[0])
    new_order = list(chain.from_iterable(step_structure))
    r = len(new_order)
    if r < n:
        new_order = new_order + list(range(r, n))
    F_out = []
    k = 0
    for i in range(0, len(step_structure)):
        tup = step_structure[i]
        reordered = reorder(F[i], new_order)
        for j in range(0, len(tup)):
            F_out.append(reordered[:k] + [F[i][x] for x in tup[:(j + 1)]] + [zero] * (n - j - 1 - k))
        k += len(tup)
    if r < n:
        for i in range(len(step_structure), m):
            F_out.append(F[i])
    return F_out


def square_forms_penalty(F, S, variants):
    F1, S1, variants1 = rearrange_penalty(F, S, variants, True)
    m = len(F1)
    if m == 0:
        return [], [], []

    # Getting the step structure
    step_structure = get_step_structure(F1)
    num_matrices = 1
    for item in step_structure:
        num_matrices = num_matrices * factorial(len(item))

    perm = map(permutations, step_structure)
    prod = product(*perm)
    stdout.write("\t" + str(num_matrices) + " matrices\n")
    stdout.flush()
    i = 0
    for order in prod:
        # if i % 100000 == 0:
        stdout.write("\tM" + str(i + 1) + "; " + str(datetime.now()) + "\n")
        stdout.flush()
        i += 1
        F_yield = squarify(F1, order)
        S_yield = squarify(S1, order)
        yield F_yield, S_yield, reorder(variants1, list(chain.from_iterable(order)))


def square_forms(F, variants):
    F1, variants1 = rearrange(F, variants, True)
    m = len(F1)
    if m == 0:
        return [], []

    # Getting the step structure
    step_structure = get_step_structure(F1)
    num_matrices = 1
    for item in step_structure:
        num_matrices = num_matrices * factorial(len(item))

    perm = map(permutations, step_structure)
    prod = product(*perm)
    stdout.write("\t" + str(num_matrices) + " matrices\n")
    stdout.flush()
    i = 0
    for order in prod:
        # if i % 100000 == 0:
        stdout.write("\tM" + str(i + 1) + "; " + str(datetime.now()) + "\n")
        stdout.flush()
        i += 1
        F_yield = squarify(F1, order)
        yield F_yield, reorder(variants1, list(chain.from_iterable(order)))


def sub_f(in_F_file, in_var_file, sub_var_file, out_F_file, out_var_file):
    f1 = open(in_F_file)
    f2 = open(in_var_file)
    f3 = open(sub_var_file)
    f4 = open(out_F_file, "w")
    f5 = open(out_var_file, "w")

    in_F = read_dm(f1)
    in_var = read_variants_list(f2)
    sub_var = read_variants_list(f3)

    in_F = list(map(list, zip(*in_F)))
    out_F = []
    out_var = []
    for i in range(0, len(in_var)):
        var = in_var[i]
        if var in sub_var:
            out_var.append(var)
            out_F.append(in_F[i])
    out_F = list(map(list, zip(*out_F)))
    out_F, out_var = rearrange(out_F, out_var, True)
    write_dm(out_F, f4)
    for var in out_var:
        f5.write(str(var) + "\n")
    f1.close()
    f2.close()
    f3.close()
    f4.close()
    f5.close()


def remap_parents(indices, parents):
    remapped = [0] * len(parents)
    for i in range(1, len(parents)):
        if parents[i] != 0:
            remapped[indices[i - 1]] = indices[parents[i] - 1]
    return remapped


def read_results(f_path):
    f = open(f_path)
    out = list(map(lambda x: list(map(int, x.split(None))), f.readlines()))
    f.close()
    return out


def read_clones(f_path):
    f = open(f_path)
    out = []
    lines = f.readlines()
    for line in lines:
        out.append(list(map(lambda x: [0] + list(map(int, x.split("-"))), line.split(";"))))
    f.close()
    return out


def write_parents(var, parents, out):
    out.write("Clone\tParent\n")
    for i in range(0, len(var)):
        out.write(str(var[i]) + "\t")
        if parents[i + 1] == 0:
            out.write("0\n")
        else:
            out.write(str(var[parents[i + 1] - 1]) + "\n")


def print_parents(var, parents):
    print("Clone\tParent")
    for i in range(0, len(var)):
        print(str(var[i]) + "\t", end=" ")
        if parents[i + 1] == 0:
            print("0")
        else:
            print(str(var[parents[i + 1] - 1]))


def white_list(variants, clones):
    out = {}
    for variant in variants:
        temp = set()
        for clone in map(set, clones):
            if variant in clone:
                temp.update(clone)
        if len(temp) == 0:
            temp = set(variants)
        temp.remove(variant)
        out[variant] = temp
    return out


def black_list(variants, clones):
    out = {}
    for variant in variants:
        out[variant] = set()
    clones_sets = list(map(set, clones))
    for i in range(0, len(clones_sets) - 1):
        clone1 = clones_sets[i]
        for j in range(i + 1, len(clones_sets)):
            clone2 = clones_sets[j]
            inter = clone1.intersection(clone2)
            if len(inter) > 0:
                sym_diff = clone1.symmetric_difference(clone2)
                for variant in sym_diff:
                    out[variant].update(inter)
    return out


def translate_clones(variants, clones):
    out = []
    for clone in clones:
        temp = [0]
        for v in clone:
            if v in variants:
                temp.append(variants.index(v) + 1)
        out.append(temp)
    return out


def sample_clones(truth, n):
    nodes = list(np.random.choice(range(1, len(truth)), replace=False, size=n))
    clones = []
    for node in nodes:
        clone = []
        temp = node
        while temp != 0:
            clone.append(temp)
            temp = truth[temp]
        clones.append(clone)
    return clones


def write_clones(clones, out_file):
    reps = []
    for clone in clones:
        clone_string = "-".join(list(map(str, clone)))
        reps.append(clone_string)
    out_file.write(";".join(reps) + "\n")


def parents_to_clones(parents):
    clones = set()
    for i in range(0, len(parents)):
        clone = set()
        clone.add(i)
        temp = parents[i]
        while temp != 0:
            # print(temp)
            clone.add(temp)
            temp = parents[temp]
        clone.add(0)
        clone = frozenset(clone)
        clones.add(clone)
    return clones


def get_clones1(parents, variants):
    variants_new = ['0'] + variants
    variants_sets = []
    for i in range(0, len(variants_new)):
        variants_sets.append(set(map(int, variants_new[i].split("_"))))
    clones = set()
    for i in range(1, len(parents)):
        clone = set()
        clone = clone.union(variants_sets[i])
        temp = parents[i]
        while temp != 0:
            clone = clone.union(variants_sets[temp])
            temp = parents[temp]
        clone.add(0)
        clone = frozenset(clone)
        clones.add(clone)
    return clones


def get_clones2(parents, removed_variants):
    variants = list(range(1, len(parents)))
    removed_variants = sorted(list(map(int, removed_variants)))
    mapping = [0] * len(parents)
    new_parents = []
    j = 0
    for i in range(0, len(parents)):
        if i not in removed_variants:
            mapping[i] = j
            parent = parents[i]
            while parent in removed_variants:
                parent = parents[parent]
            new_parents.append(mapping[parent])
            j += 1
        else:
            variants.remove(i)
    variants = list(map(str, variants))
    return get_clones1(new_parents, variants)
