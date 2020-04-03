import numpy as np
from myutils import *
from math import *
from predict import *
from time import *
import sys


def split_one(n):
    out = list(map(Decimal, np.random.rand(n)))
    s = sum(out)
    for i in range(n):
        out[i] = out[i]/s
    return out


def split_one_distributed(prev_row, n):
    out = list(map(Decimal, np.random.dirichlet([i * 10 for i in prev_row[:n]])))
    for i in range(0, len(out)):
        if out[i] == zero:
            out[i] = Decimal("0.0001")
    out.append(Decimal(np.random.uniform(0.1, 0.2)))
    s = sum(out)
    for i in range(n+1):
        out[i] = out[i] / s
    return out


def simple_simulate_square(n):
    parents = [0] * n
    C = [[0] * n]
    C[0][0] = 1
    indices = list(range(n))
    likelihood = one
    for i in range(1, n):
        parent = np.random.choice(indices, p=C[i-1])
        parents[i] = parent
        likelihood = likelihood * C[i-1][parent]
        # new_row = split_one(i+1) + ([zero] * (n - i - 1))
        new_row = split_one_distributed(C[i-1], i) + ([zero] * (n - i - 1))
        C.append(new_row)
    return C, parents, likelihood


def simple_simulator(params):
    [n, p] = params
    while True:
        C, parents, likelihood = simple_simulate_square(n)
        n_time_points = np.random.binomial(n, p)
        if n_time_points < 2:
            continue
        chosen_time_points = list(np.random.choice(range(1, n-1), replace=False, size=n_time_points-1))
        chosen_time_points.append(n-1)
        C_chosen = []
        for i in sorted(chosen_time_points):
            C_chosen.append(C[i])
        T = list_to_tree(parents)
        F = get_f(C_chosen, T)
        for row in F:
            del row[0]
        if not all([v == 0 for v in parents]):
            yield F, parents, likelihood


def write_simulations(simulator, params, N, sim_file):
    sim_file_f = open(sim_file, "w")
    sim = simulator(params)
    for i in range(N):
        F, parents, likelihood = sim.__next__()
        sim_file_f.write("# " + str(i + 1) + " " + str(len(F)) + "\n")
        sim_file_f.write("Parents: " + " ".join(map(str, parents)) + "\n")
        sim_file_f.write("Likelihood: " + str(likelihood) + "\n")
        write_dm(F, sim_file_f)
    sim_file_f.close()


def read_simulations(sim_file):
    sim_file_f = open(sim_file)
    F_list = []
    parents_list = []
    likelihood_list = []
    R_list = []
    lines = sim_file_f.readlines()
    sim_file_f.close()
    if lines[0][0] != '#':
        sys.exit("Invalid input file!")
    i = 0
    while i < len(lines):
        line = lines[i]
        words = line.split(None)
        num_times = int(words[2])
        parents = list(map(int, lines[i + 1][9:].split(None)))
        likelihood = Decimal(lines[i + 2][12:])
        F = []
        for j in range(i + 3, i + 3 + num_times):
            F.append(list(map(Decimal, lines[j].split(None))))
        R = []
        if len(words) > 3:
            for j in range(i + 3 + num_times, i + 3 + num_times + num_times):
                R.append(list(map(Decimal, lines[j].split(None))))
        F_list.append(F)
        parents_list.append(parents)
        likelihood_list.append(likelihood)
        R_list.append(R)
        if len(words) > 3:
            i = i + 3 + num_times + num_times
        else:
            i = i + 3 + num_times
    return F_list, parents_list, likelihood_list, R_list


def read_output(output_file):
    out = []
    output_file_f = open(output_file)
    lines = output_file_f.readlines()
    output_file_f.close()
    if lines[0][0] != '#':
        sys.exit("Invalid input file!")
    i = 0
    while i < len(lines):
        parents = list(map(int, lines[i + 1][9:].split(None)))
        variants = list(map(str, lines[i + 2][10:].split(None)))
        score = Decimal(lines[i + 3][7:])
        running_time = Decimal(lines[i + 4][14:])
        num_tp = int(lines[i + 5][17:])
        removed_variants = list(map(str, lines[i + 6][9:].split(None)))
        out.append((parents, variants, score, running_time, num_tp, removed_variants))
        i = i + 7
    return out


def evaluate_batch(sim_file, algos):
    F_list, parents_list, likelihood_list, R_list = read_simulations(sim_file)

    match_percentages_dict = {}
    avg_number_of_true_clones_dict = {}
    log_running_times_dict = {}
    avg_number_of_time_points_dict = {}
    for j in algos:
        match_percentages = []
        number_of_true_clones = []
        running_times = []
        number_of_time_points = []
        predicted_output = read_output(sim_file + "." + str(j))
        for i in range(0, len(parents_list)):
            true_clones = get_clones2(parents_list[i], predicted_output[i][5])
            if predicted_output[i][0]:
                predicted_clones = get_clones1(predicted_output[i][0], predicted_output[i][1])
                match_count = len(true_clones.intersection(predicted_clones))
                match_percentage = (match_count + 0.0) / len(true_clones)
            else:
                match_percentage = 0.0
            match_percentages.append(match_percentage)
            number_of_true_clones.append(len(true_clones))
            running_times.append(predicted_output[i][3])
            number_of_time_points.append(predicted_output[i][4])
        match_percentages_dict[j] = match_percentages
        avg_number_of_true_clones_dict[j] = (sum(number_of_true_clones) + 0.0) / len(number_of_true_clones)
        running_times = list(map(lambda x: x if x != 0 else 0.0001, running_times))
        log_running_times = list(map(np.log10, running_times))
        log_running_times_dict[j] = log_running_times
        avg_number_of_time_points_dict[j] = (sum(number_of_time_points) + 0.0) / len(number_of_time_points)

    o1 = open(sim_file + ".recall", "w")
    for i in range(0, len(parents_list)):
        for j in algos:
            o1.write("{0:.5f}".format(match_percentages_dict[j][i]) + "\t")
        o1.write("\n")
    o1.close()

    o2 = open(sim_file + ".info", "w")
    for i in algos:
        o2.write(str(avg_number_of_true_clones_dict[i]) + "\t")
    o2.write("\n")
    for i in algos:
        o2.write(str(avg_number_of_time_points_dict[i]) + "\t")
    o2.write("\n")
    o2.close()

    o3 = open(sim_file + ".ltimes", "w")
    for i in range(0, len(parents_list)):
        for j in algos:
            o3.write("{0:.5f}".format(log_running_times_dict[j][i]) + "\t")
        o3.write("\n")
    o3.close()


def predict_batch(sim_file, algo):
    stdout.write("\n" + str(datetime.now()) + "\n")
    stdout.write("Algo " + str(algo) + "\n")
    stdout.flush()
    out_file_f = open(sim_file + "." + str(algo), "w")
    F_list, parents_list, likelihood_list, R_list = read_simulations(sim_file)
    for i in range(0, len(F_list)):
        stdout.write("Sim " + str(i + 1) + "; " + str(datetime.now()) + "\n")
        stdout.flush()
        F = F_list[i]
        R = R_list[i]
        var = list(range(1, len(F[0]) + 1))
        parents, score, variants, removed_variants, num_times, running_time, _ = predict(F, var, algo, R)
        out_file_f.write("# " + str(i + 1) + "\n")
        out_file_f.write("Parents: " + " ".join(map(str, parents)) + "\n")
        out_file_f.write("Variants: " + " ".join(map(str, variants)) + "\n")
        out_file_f.write("Score: " + str(score) + "\n")
        out_file_f.write("Running Time: " + str(running_time) + "\n")
        out_file_f.write("Num time-points: " + str(num_times) + "\n")
        out_file_f.write("Removed: " + " ".join(map(str, removed_variants)) + "\n")
    out_file_f.close()


def introduce_noise(F):
    F_out = []
    R_out = []
    rd_mean = 170
    rd_sd = 40
    for i in range(0, len(F)):
        row = []
        R_row = []
        for j in range(0, len(F[i])):
            if F[i][j] == 0:
                row.append(F[i][j])
                R_row.append(zero)
            else:
                random_rd = Decimal(abs(np.random.normal(rd_mean, rd_sd)))
                sd = Decimal(sqrt((F[i][j] * (1 - F[i][j])) / random_rd))
                new_val = Decimal(np.random.normal(F[i][j], sd))
                if new_val < zero:
                    new_val = Decimal("0.00001")
                elif new_val > one:
                    new_val = one
                row.append(new_val)
                R_row.append(random_rd)
        F_out.append(row)
        R_out.append(R_row)
    return F_out, R_out


def add_noise_to_simulation(sim_file):
    out_sim_file = open(sim_file + ".noisy", "w")
    F_list, parents_list, likelihood_list, _ = read_simulations(sim_file)
    for i in range(0, len(F_list)):
        F_noisy, R = introduce_noise(F_list[i])
        out_sim_file.write("# " + str(i + 1) + " " + str(len(F_noisy)) + " noisy\n")
        out_sim_file.write("Parents: " + " ".join(map(str, parents_list[i])) + "\n")
        out_sim_file.write("Likelihood: " + str(likelihood_list[i]) + "\n")
        write_dm(F_noisy, out_sim_file)
        write_dm(R, out_sim_file)
    out_sim_file.close()



