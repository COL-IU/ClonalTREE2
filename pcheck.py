from myutils import *
import copy


class to_return:  # class consiting a number and a list
    def __init__(self, maxprob, bestpar):
        self.maxprob = maxprob
        self.bestpar = bestpar

    def __str__(self):
        return "The probability: " + str(self.maxprob) + " The parents: " + str(self.bestpar)


def general_readin(inputfile):  # read in from a file
    mymatrix = []
    with open(inputfile, 'r') as f:  # read in the matrix
        for line in f:
            if len(line) > 0:
                mymatrix.append(line.split())
    for i in range(len(mymatrix)):
        for j in range(len(mymatrix[i])):
            mymatrix[i][j] = Decimal(mymatrix[i][j])  # convert the elements into numbers
    return mymatrix


def line_pre_check(vec, par):  # precheck in the given order(simple case) returns the parent options
    for i in range(len(par)):
        if (par[i] != -1):
            vec[par[i]] = vec[par[i]] - vec[i + 1]
    modified = True
    while (modified):
        parentoptions = list()
        modified = False
        for i in range(1, len(vec)):
            if (par[i - 1] != -1):
                parentoptions.append([par[i - 1]])
            if (par[i - 1] == -1):
                options = 0
                parent = -1
                paropt = list()
                for j in range(i):
                    if (vec[j] >= vec[i]):
                        options = options + 1
                        parent = j
                        paropt.append(j)
                parentoptions.append(paropt)
                if (options == 0):
                    return -1
                if (options == 1):
                    par[i - 1] = parent
                    vec[parent] = vec[parent] - vec[i]
                    modified = True
        return parentoptions


def matrix_pre_check(matrix, par):  # call line_pre_check for each row_index
    optionalparents = []
    matrixcopy = copy.deepcopy(matrix)
    for i in matrixcopy:
        lineoptions = line_pre_check(i, par)
        optionalparents.append(lineoptions)
    return optionalparents


def full_pre_check(matrix, par):  # calling matrix_pre_check until there is no change
    parcopy = copy.deepcopy(par)
    newpar = []
    listoflists = matrix_pre_check(matrix, par)
    for i in range(len(par)):  # merging the different conditions for a given node
        intersect = set(listoflists[0][i])
        for j in listoflists:
            if j != -1:
                intersect = intersect.intersection(j[i])
        if (len(intersect) == 1):
            newpar.append(list(intersect)[0])  # label the parent when only one option
        else:
            newpar.append(-1)
    if (newpar.count(-1) != parcopy.count(-1)):  # if new information call again
        newpar = full_pre_check(matrix, newpar)
    return newpar


def count_prob(vec,
               par):  # with given directed tree returns the number we want to maximise and whether the tree is appropiate
    if ((len(vec) - 1) != len(par)):  # returns -1 if not appropriate
        return ("Not appropiate length")
    prob = 1
    good = True
    for i in range(len(par)):
        x = vec[par[i]]
        for j in range(i):
            if (par[j] == par[i]):
                x = x - vec[j + 1]
        prob = prob * x
        good = (good and (x - vec[i + 1]) >= 0)
    if (not good):
        x = to_return(-1, [])
        return x
    x = to_return(prob, par)
    return x


def general_count_prob_good(matrix, par):  # Right way: we use the former row_index only for the first newborn
    for i in matrix:
        if (count_prob(i, par).maxprob == -1):
            return to_return(-1, [])
    newclones = []
    for i in matrix:
        newc = 0
        for j in i:
            if (j != 0):
                newc = newc + 1
        newclones.append(newc)  # how many clones in each row_index
    for i in range(len(newclones) - 1, 0, -1):
        newclones[i] = newclones[i] - newclones[i - 1]  # how many clones new in a given row_index
    k = 0
    l = 0
    lowertriangle = []
    for i in range(len(newclones)):
        for j in range(newclones[i]):
            lowertriangle.append(matrix[k][0:l + 1])
            l = l + 1
        k = k + 1
    for i in range(len(lowertriangle)):
        for j in range(len(lowertriangle)):
            if (j > i):
                lowertriangle[i].append(0)
    prob = 1
    for i in range(1, len(lowertriangle)):
        actprob = lowertriangle[i - 1][par[i - 1]]
        for j in range(i - 1):
            if (par[j] == par[i - 1]):
                actprob = actprob - lowertriangle[i - 1][j + 1]
        prob = prob * actprob
    y = to_return(prob, par)
    return y


def general_recursive(origmatrix, vec, par, indsort):  # general recursive algorithm used in the general case
    maxiprob = to_return(-2, [])
    for i in range(1, len(indsort)):
        if (par[indsort[i] - 1] != -1):
            continue
        else:
            maxiprob = to_return(-1, [])
            for j in range(i):
                if ((indsort[j] < indsort[i]) and (vec[indsort[j]] >= vec[indsort[i]])):
                    par2 = list(par)
                    vec2 = list(vec)
                    par2[indsort[i] - 1] = indsort[j]
                    vec2[indsort[j]] = vec2[indsort[j]] - vec2[indsort[i]]
                    x = general_recursive(origmatrix, vec2, par2, indsort)
                    if (x.maxprob > maxiprob.maxprob):
                        maxiprob.maxprob = x.maxprob
                        maxiprob.bestpar = x.bestpar
            break
    if (maxiprob.maxprob == -2):  # not changed to minus 1, it means the parent vector is fully filled
        maxiprob = general_count_prob_good(origmatrix, par)
    return maxiprob


# tempmatrix = general_readin('input5.txt')
# mymatrix = add_founder(tempmatrix)
# par = [-1] * (len(mymatrix) - 1)
#
# knownparents = full_pre_check(mymatrix, par)
# lastrow = mymatrix[len(mymatrix) - 1]
# lastrowcopy = copy.deepcopy(lastrow)
# indsort = sorted(range(len(lastrowcopy)), key=lambda k: lastrowcopy[k], reverse=True)
# for i in range(len(knownparents)):
#     if (knownparents[i] != -1):
#         lastrowcopy[knownparents[i]] = lastrowcopy[knownparents[i]] - lastrowcopy[i + 1]
#
# print(knownparents)
# print(general_recursive(mymatrix, lastrowcopy, knownparents, indsort))
