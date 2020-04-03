from myutils import *
import numpy as np
from time import *
from pcheck import *
import sys


def filter_data(F, R, variants):
    m = len(F)
    n = len(F[0])
    F1 = [[zero for _ in range(len(F[0]))] for _ in range(len(F))]
    R1 = [[zero for _ in range(len(R[0]))] for _ in range(len(R))]
    for i in range(0, m):
        for j in range(0, n):
            if F[i][j] >= 0.05:
                F1[i][j] = F[i][j]
                R1[i][j] = R[i][j]
            else:
                F1[i][j] = zero
                R1[i][j] = zero
    F1 = list(map(list, zip(*F1)))
    R1 = list(map(list, zip(*R1)))
    chosen_variants = []
    removed_variants = []
    F2 = []
    R2 = []
    for i in range(0, n):
        seen_times = []
        for j in range(0, m):
            if F1[i][j] > zero:
                seen_times.append(j)
        if len(seen_times) >= 2 and seen_times == list(range(min(seen_times), max(seen_times) + 1)):
            chosen_variants.append(variants[i])
            F2.append(F1[i])
            R2.append(R1[i])
        else:
            removed_variants.append(variants[i])
    F2 = list(map(list, zip(*F2)))
    R2 = list(map(list, zip(*R2)))
    if chosen_variants:
        F3, R3, rearranged_variants, _ = rearrange_penalty(F2, R2, chosen_variants)
        return F3, R3, rearranged_variants, removed_variants
    else:
        return F2, R2, chosen_variants, removed_variants


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


def cluster_data(F, R, variants, threshold):
    F_t = list(map(list, zip(*F)))
    R_t = list(map(list, zip(*R)))
    str_variants = list(map(str, variants))
    clusters = [(str_variants[0], F_t[0], R_t[0])]
    i = 1
    while i < len(F_t):
        new_cluster = True
        for j in range(0, len(clusters)):
            (header, column, r_column) = clusters[j]
            if diff_below_cutoff(F_t[i], column, threshold):
                clusters[j] = (header + "_" + str_variants[i], average_column(F_t[i], column), average_column(R_t[i], r_column))
                new_cluster = False
                break
        if new_cluster:
            clusters.append((str_variants[i], F_t[i], R_t[i]))
        i = i + 1
    cluster_variants = []
    cluster_F = []
    cluster_R = []
    for (cluster_variant, column, r_column) in clusters:
        cluster_variants.append(cluster_variant)
        cluster_F.append(column)
        cluster_R.append(r_column)
    cluster_F = list(map(list, zip(*cluster_F)))
    cluster_R = list(map(list, zip(*cluster_R)))
    return cluster_F, cluster_R, cluster_variants


def best_scoring_ancestor(F, S, order, parents, clone_index, prev_clone_arrival_time, curr_clone_arrival_time):
    num_clones = len(order)
    scores = [zero] * len(order)
    row_index = prev_clone_arrival_time

    for i in range(0, clone_index):
        clone = order[i]
        scores[clone] = F[row_index][clone]

    children = parents_to_children(parents, num_clones)
    new_scores = [zero] * len(order)
    for i in range(0, clone_index):
        clone = order[i]
        children_sum = zero
        cur_children = children[clone]
        for child in cur_children:
            children_sum = children_sum + F[row_index][child]
        new_scores[clone] = scores[clone] - children_sum

    penalized_scores = [zero] * len(order)
    for i in range(0, clone_index):
        clone = order[i]
        cur_children = children[clone]
        cur_children.append(order[clone_index])
        penalty = zero
        for t in range(curr_clone_arrival_time, len(F)):
            children_sum = zero
            cs_sds = zero
            for child in cur_children:
                children_sum = children_sum + F[t][child]
                cs_sds = cs_sds + (S[t][child] ** 2)
            temp = F[t][clone] - children_sum
            if temp < 0:
                denominator = Decimal(sqrt((S[t][clone] ** 2) + cs_sds))
                if denominator > zero:
                    penalty = penalty + (temp / denominator)
                    # penalty = penalty + temp
        penalized_scores[clone] = new_scores[clone] + penalty

    temp = []
    for i in range(0, clone_index):
        temp.append((order[i], penalized_scores[order[i]]))
    ancestor = max(temp, key=lambda i: i[1])[0]

    return ancestor, new_scores[ancestor], penalized_scores[ancestor] - new_scores[ancestor]


def tree_score(F, S, order, parents, arrival_times, upto):
    num_clones = len(order)
    children = parents_to_children(parents, num_clones)
    final_score = one
    final_penalty = zero
    for clone_index in range(1, upto + 1):
        prev_clone = order[clone_index - 1]
        curr_clone = order[clone_index]
        parent = parents[curr_clone]
        prev_clone_arrival_time = arrival_times[prev_clone]
        curr_clone_arrival_time = arrival_times[curr_clone]

        row_index = prev_clone_arrival_time
        score = F[row_index][parent]

        children_sum = zero
        cur_children = children[parent]
        for child in cur_children:
            children_sum = children_sum + F[row_index][child]
        new_score = score - children_sum

        penalty = zero
        for t in range(curr_clone_arrival_time, len(F)):
            children_sum = zero
            cs_sds = zero
            for child in cur_children:
                children_sum = children_sum + F[t][child]
                cs_sds = cs_sds + (S[t][child] ** 2)
            temp = F[t][parent] - children_sum
            if temp < 0:
                denominator = Decimal(sqrt((S[t][parent] ** 2) + cs_sds))
                if denominator > zero:
                    # print(str(temp) + "\t" + str(denominator) + "\t" + str(temp / denominator))
                    penalty = penalty + (temp / denominator)
                    # penalty = penalty + temp
        final_penalty = final_penalty + penalty
        final_score = final_score * new_score

    return final_score, final_penalty


def new_algorithm(F, R):
    # Empty check
    if len(F) == 0:
        return [], -1  # parents_vector, score

    S = rd_to_sd(F, R)

    # Add founder
    my_F = add_founder(F)
    my_S = add_founder_penalty(S)
    num_clones = len(my_F[0])

    # Get the step structure of the step matrix
    steps, arrival_times = get_step_structure(my_F)

    parents = {}
    for clone in range(1, num_clones):
        parents[clone] = -1
    order = list(range(0, num_clones))
    final_score = one
    final_penalty = zero
    clone_counter = 1
    for step_counter in range(1, len(steps)):
        step = steps[step_counter]
        if len(step) == 1:
            clone = order[clone_counter]
            parent, score, penalty = best_scoring_ancestor(my_F, my_S, order, parents, clone_counter,
                                                           arrival_times[order[clone_counter - 1]],
                                                           arrival_times[order[clone_counter]])
            parents[clone] = parent
            final_score = final_score * score
            final_penalty = final_penalty + penalty
            clone_counter = clone_counter + 1
        else:
            G, vertices = get_partial_order(my_F, [step[0], step[-1]])
            ts = topological_sort(G, vertices)
            new_order = order[:step[0]] + ts + order[step[-1] + 1:]
            for i in range(0, len(ts)):
                clone = new_order[clone_counter]
                parent, score, penalty = best_scoring_ancestor(my_F, my_S, new_order, parents, clone_counter,
                                                               arrival_times[order[clone_counter - 1]],
                                                               arrival_times[order[clone_counter]])
                max_parents = parents.copy()
                max_parents[clone] = parent
                # print(scores[parent])
                max_score = final_score * score
                max_penalty = final_penalty + penalty
                max_order = new_order[:]
                for j in range(0, i):
                    new_ts = ts[:]
                    new_ts.remove(ts[i])
                    new_ts.insert(j, ts[i])
                    common_ancestor = parents[new_ts[j + 1]]
                    temp_parents = parents.copy()
                    temp_parents[clone] = common_ancestor
                    temp_order = order[:step[0]] + new_ts + order[step[-1] + 1:]
                    temp_score, temp_penalty = tree_score(my_F, my_S, temp_order, temp_parents, arrival_times,
                                                          clone_counter)
                    if (temp_score + temp_penalty) > (max_score + max_penalty):
                        max_parents = temp_parents.copy()
                        max_score = temp_score
                        max_penalty = temp_penalty
                        max_order = temp_order
                parents = max_parents.copy()
                final_score = max_score
                final_penalty = max_penalty
                new_order = max_order[:]
                clone_counter = clone_counter + 1
            order = new_order[:]
    parents_vector = [0]
    for i in range(1, num_clones):
        parents_vector.append(parents[i])

    return parents_vector, final_score + final_penalty


def valid_parent_value(parents, F, fail_threshold):
    # F may not be square
    m = len(F)
    vp = True
    cur_parent = parents[-1]
    children = get_children(parents, cur_parent)
    for t in range(cur_parent, m):
        children_sum = zero
        for child in children:
            children_sum += F[t][child]
        if F[t][cur_parent] - children_sum < fail_threshold:
            vp = False
            return vp
    return vp


def valid_parent_order(order_validity, parents):
    cur_variant = len(parents)-1
    cur_path = deepcopy(order_validity[cur_variant])
    temp = parents[cur_variant]
    while temp != 0:
        cur_path.discard(temp)
        temp = parents[temp]
    cur_path.discard(0)
    return len(cur_path) == 0


def c_row(f_row, parents):
    m = len(f_row)
    C = []
    for parent in range(0, m):
        children = get_children(parents, parent)
        children_sum = zero
        for child in children:
            children_sum += f_row[child]
        C.append(f_row[parent] - children_sum)
    return C


def smart_predict_original(F, m, fail_threshold, clones=[], S=[]):
    # Assumes square F
    if len(F) == 0:
        return [], [], [], -1

    # Add founder
    my_F = add_founder(F)
    # m = len(my_F)
    choices_stack = [[0]]
    chosen_parents = [0]

    valid_parents = white_list(list(range(m)), clones)
    order_validity = black_list(list(range(m)), clones)

    success = False
    while choices_stack:
        chosen_parents = choices_stack.pop()
        i = len(chosen_parents) - 1
        if not valid_parent_value(chosen_parents, my_F, fail_threshold):
            continue
        if clones:
            if not valid_parent_order(order_validity, chosen_parents):
                continue
        if i == (m - 1):
            success = True
            break
        C_row = c_row(my_F[i], chosen_parents)
        next_choices = list((np.array(C_row)).argsort())
        for next_choice in next_choices:
            if C_row[next_choice] > fail_threshold and next_choice <= i and next_choice in valid_parents[i+1]:
                temp = chosen_parents[:]
                temp.append(next_choice)
                choices_stack.append(temp)
    if not success:
        return [], my_F, [], -1

    C = []
    for i in range(0, m):
        C.append(c_row(my_F[i], chosen_parents))

    p = one
    for i in range(1, m):
        p = p * C[i-1][chosen_parents[i]]
    return C, my_F, chosen_parents, p


def variant_penalty1(parents, F, S):
    # F may not be square
    m = len(F)
    pen = 0
    cur_parent = parents[-1]
    children = get_children(parents, cur_parent)
    for t in range(cur_parent, m):
        children_sum = zero
        cs_sds = zero
        for child in children:
            children_sum += F[t][child]
            cs_sds += (S[t][child] ** 2)
        temp = F[t][cur_parent] - children_sum
        if temp < 0:
            denominator = Decimal(sqrt((S[t][cur_parent] ** 2) + cs_sds))
            if denominator > zero:
                pen += temp / denominator
    return pen


def scores(F, parents, penalty, S=[]):
    i = len(parents) - 1
    C_row = c_row(F[i], parents)
    for j in range(0, len(C_row)):
        if S:
            pen = penalty(parents + [j], F, S)
            C_row[j] = C_row[j] + pen
    return C_row


def smart_predict_penalty(F, m, fail_threshold, clones=[], S=[]):
    # Assumes square F
    if len(F) == 0:
        return [], [], [], -1

    # Add founder
    my_F = add_founder(F)
    my_S = []
    if S:
        my_S = add_founder_penalty(S)
    # m = len(my_F)
    choices_stack = [[0]]
    chosen_parents = [0]

    valid_parents = white_list(list(range(m)), clones)
    order_validity = black_list(list(range(m)), clones)

    success = False
    while choices_stack:
        # print(choices_stack)
        chosen_parents = choices_stack.pop()
        i = len(chosen_parents) - 1
        if not S:
            if not valid_parent_value(chosen_parents, my_F, fail_threshold):
                continue
        if clones:
            if not valid_parent_order(order_validity, chosen_parents):
                continue
        if i == (m - 1):
            success = True
            break
        sc = scores(my_F, chosen_parents, variant_penalty1, my_S)
        next_choices = list((np.array(sc)).argsort())
        for next_choice in next_choices:
            if next_choice <= i and next_choice in valid_parents[i + 1]:
                temp = chosen_parents[:]
                temp.append(next_choice)
                choices_stack.append(temp)
    if not success:
        return [], my_F, [], -1

    C = []
    for i in range(0, m):
        C.append(c_row(my_F[i], chosen_parents))

    p = one
    for i in range(1, m):
        p = p * C[i-1][chosen_parents[i]]
    return C, my_F, chosen_parents, p


def old_algorithm(F, R, smart_predict_algo):
    # Empty check
    if len(F) == 0:
        return [], -1  # parents_vector, score

    S = rd_to_sd(F, R)

    num_clones = len(F[0])

    steps, arrival_times = get_step_structure(F)
    choices_stack = [(0, [0], [])]
    (cur_score, cur_T, cur_P) = (0, [0], [])

    success = False
    while choices_stack:
        (cur_score, cur_T, cur_P) = choices_stack.pop()
        cur_m = 1
        for item in cur_P:
            cur_m += len(item)
        i = len(cur_P)
        if i == len(steps):
            success = True
            break
        cur_step = steps[i]
        temp = []
        j = 0
        # topo_permutations, F_mod = topological_valid_permutations(F_iter, cur_step)
        # F_iter = deepcopy(F_mod)
        # for permutation in topo_permutations:
        for permutation in permutations(cur_step):
            perm_ss = deepcopy(cur_P)
            perm_ss.append(list(permutation))
            perm_F = squarify(F, perm_ss)
            perm_S = squarify(S, perm_ss)
            _, _, perm_T, perm_score = smart_predict_algo(perm_F, cur_m + len(permutation), zero, clones=[], S=perm_S)
            if perm_score != -1:
                temp.append((perm_score, perm_T, perm_ss))
        if temp:
            temp_sorted = sorted(temp, key=lambda x: x[0])
            choices_stack = choices_stack + temp_sorted

    if success:
        variants_in = list(range(1, num_clones + 1))
        variants_out = reorder(variants_in, chain.from_iterable(cur_P))
        parents_out = remap_parents(variants_out, cur_T)
        return parents_out, cur_score
    else:
        return [], -1


def predict(F, variants, algo, R=[], filter=True):
    start = process_time()
    if len(R) == 0:
        R = [[Decimal(sys.maxsize) for _ in range(len(F[0]))] for _ in range(len(F))]
    if filter:
        F, R, variants, removed_variants = filter_data(F, R, variants)
    else:
        removed_variants = []
    F, R, _, removed_time_points = remove_redundant_time_points(F, R)

    if variants:
        # F, R, variants = cluster_data(F, R, variants, 0.05)
        # F, R, _, _ = remove_redundant_time_points(F, R)
        parents = []
        score = zero
        if algo == 0:
            parents, score = old_algorithm(F, R, smart_predict_original)
        elif algo == 1:
            parents, score = new_algorithm(F, R)
        elif algo == 2:
            parents, score = old_algorithm(F, R, smart_predict_penalty)
        else:
            exit("Invalid parameter for algo.")
    else:
        parents = []
        score = zero
    end = process_time()

    return parents, score, variants, removed_variants, len(F), end-start, removed_time_points
