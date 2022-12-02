import helper as my
import motifSearch as I
from functools import lru_cache
import itertools
import math
import random
import copy
import networkx as nx
import scipy
import pickle

# =============================================================================
# SECTION 1
# =============================================================================
def binomial_coeff(m,n):
    silnia = scipy.math.factorial
    return (int(silnia(m+n)/(float(silnia(m)*silnia(n)))))


def greedy_change(money):
    coins = []
    denomination = [100, 50, 25, 10, 5, 1]
    while money > 0:
        coin = [x for x in denomination if x <= money][0]
        coins.append(coin)
        money -= coin
    return coins


@lru_cache(maxsize=128, typed=False) # cache, no need to recalculate
def recursive_change(money, coins):  # NB missing def dpchange()
    if money == 0:
        return 0
    min_num_coins = float('inf')
    for coin in coins:
        if money >= coin:
            num_coins = recursive_change(money-coin, coins)
            if num_coins + 1  < min_num_coins:
                min_num_coins = num_coins + 1
    return min_num_coins


def manhattan_tourist(n, m, down, right):
    s = [([0]*(m+1)) for i in range(n+1)]
    for i in range(n+1):
        for j in range(m+1):
            if i == 0 and j == 0:
                continue
            if i == 0:
                s[0][j] = s[0][j-1] + right[0][j-1]
            elif j == 0:
                s[i][0] = s[i-1][0] + down[i-1][0]
            else:
                s[i][j] = max(s[i-1][j] + down[i-1][j],
                              s[i][j-1] + right[i][j-1])
    return s[n][m]


def lcs(xs, ys):
    if xs and ys:
        *xb, xe = xs
        *yb, ye = ys
        if xe == ye:
            return lcs(xb, yb) + [xe]
        else:
            return max(lcs(xs, yb), lcs(xb, ys), key=len)
    else:
        return []


def lcs_backtrack(v, w):
    n = len(v)
    m = len(w)
    s = [([0]*(m+1)) for i in range(n+1)]
    backtrack = [([' ']*(m+1)) for i in range(n+1)]
    for i in range(n+1):
        for j in range(m+1):
            if i == 0 or j == 0:
                continue
            else:
                match = 0
                if v[i-1] == w[j-1]:
                    match = 1
                s[i][j] = max(s[i-1][j], s[i][j-1], s[i-1][j-1] + match)
                if s[i][j] == s[i-1][j]:
                    backtrack[i][j] = '↓'
                elif s[i][j] == s[i][j-1]:
                    backtrack[i][j] = '→'
                elif s[i][j] == s[i-1][j-1] + match:
                    backtrack[i][j] = '↘'

    
    def lcs_output(backtrack, v, i, j):
        if i == 0 or j == 0:
            return ''
        if backtrack[i][j] == '↓':
            return lcs_output(backtrack, v, i-1, j)
        elif backtrack[i][j] == '→':
            return lcs_output(backtrack, v, i, j-1)
        else:
            return lcs_output(backtrack, v, i-1, j-1) + v[i-1]
        
    return lcs_output(backtrack, v, n, m)


def dag_formatting(data):
    us, vs, ws = [], [], []
    for line in data:
        temp = line.split('->')
        us.append((temp[0]))
        vs.append((temp[1].split(':')[0]))
        ws.append(int(temp[1].split(':')[1]))
        
    G = nx.DiGraph()
    
    for u, v, w in zip(us, vs, ws):
        G.add_edge(u, v, weight=w)

    return G


def dag_long_path(G):
    return nx.dag_longest_path_length(G), nx.dag_longest_path(G)


def aligner(backtrack, a, b, n, m, mark):
    align_a, align_b = '', ''
    while n+1 != 0 and m+1 != 0:
        if backtrack[n][m] == mark[0]:
            align_a += a[n-1]
            align_b += b[m-1]
            n -= 1
            m -= 1
        elif backtrack[n][m] == mark[1]:
            align_a += a[n-1]
            align_b += '-'
            n -= 1
        elif backtrack[n][m] == mark[2]:
            align_a += '-'
            align_b += b[m-1]
            m -= 1
        else:
            break    
    return align_a[::-1], align_b[::-1]


def prot_global_align(a, b, gap=5):
    with open('data\sm\BLOSUM62.pkl', 'rb') as input_file:
        score = pickle.load(input_file)
    
    n, m = len(a), len(b)
    s = [([0]*(m+1)) for i in range(n+1)]
    backtrack = [([' ']*(m+1)) for i in range(n+1)]
    mark = ['↘', '↓', '→']

    for i in range(n+1):
        for j in range(m+1):
            if i == 0 and j == 0:
                continue
            if i == 0:
                backtrack[i][j] = mark[2]
                s[i][j] = s[0][j-1] - gap
            elif j == 0:
                backtrack[i][j] = mark[1]
                s[i][j] = s[i-1][0] - gap
            else:
                test = [s[i-1][j-1] + score[a[i-1]][b[j-1]],
                        s[i-1][j] - gap,
                        s[i][j-1] - gap]
                s[i][j] = max(test)
                
                index = test.index(max(test))
                backtrack[i][j] = mark[index]
    
    max_score = s[n][m]
    
    align_a, align_b = aligner(backtrack, a, b, n, m, mark)
    
    return max_score, align_a, align_b


def prot_local_align(a, b, gap=5):
    with open('data\sm\Pam250.pkl', 'rb') as input_file:
        score = pickle.load(input_file)
    
    n, m = len(a), len(b)
    s = [([0]*(m+1)) for i in range(n+1)]
    
    backtrack = [([' ']*(m+1)) for i in range(n+1)]
    mark = ['↘', '↓', '→', '*']

    for i in range(n+1):
        for j in range(m+1):
            if i == 0 and j == 0:
                continue
            if i == 0:
                backtrack[i][j] = mark[2]
                s[i][j] = 0
            elif j == 0:
                backtrack[i][j] = mark[1]
                s[i][j] = 0
            else:
                test = [s[i-1][j-1] + score[a[i-1]][b[j-1]],
                        s[i-1][j] - gap,
                        s[i][j-1] - gap,
                        0]
                s[i][j] = max(test)
                
                index = test.index(max(test))
                backtrack[i][j] = mark[index]
    
    max_score, max_index = max((x, (i, j))
                               for i, row in enumerate(s)
                               for j, x in enumerate(row))
    
    n, m = max_index[0], max_index[1]
    
    align_a, align_b = '', ''
    while backtrack[n][m] != mark[3]:
        if backtrack[n][m] == mark[0]:
            align_a += a[n-1]
            align_b += b[m-1]
            n -= 1
            m -= 1
        elif backtrack[n][m] == mark[1]:
            align_a += a[n-1]
            align_b += '-'
            n -= 1
        elif backtrack[n][m] == mark[2]:
            align_a += '-'
            align_b += b[m-1]
            m -= 1
    
    return max_score, align_a[::-1], align_b[::-1]


def edit_distance(a, b, gap=1):
    n, m = len(a), len(b)
    s = [([0]*(m+1)) for i in range(n+1)]
    backtrack = [([' ']*(m+1)) for i in range(n+1)]
    mark = ['↘', '↓', '→']

    for i in range(n+1):
        for j in range(m+1):
            if i == 0 and j == 0:
                continue
            if i == 0:
                backtrack[i][j] = mark[2]
                s[i][j] = s[0][j-1] - gap
            elif j == 0:
                backtrack[i][j] = mark[1]
                s[i][j] = s[i-1][0] - gap
            else:
                match = 0
                if a[i-1] != b[j-1]:
                    match = gap
                test = [s[i-1][j-1] - match,
                        s[i-1][j] - gap,
                        s[i][j-1] - gap]
                s[i][j] = max(test)
                
                index = test.index(max(test))
                backtrack[i][j] = mark[index]
    
    align_a, align_b = aligner(backtrack, a, b, n, m, mark)
    
    return I.hamming_distance(align_a, align_b)


def fitting_align(a, b, gap=1):
    n, m = len(a), len(b)
    s = [([0]*(m+1)) for i in range(n+1)]
    backtrack = [([' ']*(m+1)) for i in range(n+1)]
    mark = ['↘', '↓', '→']

    for i in range(n+1):
        for j in range(m+1):
            if i == 0 and j == 0:
                continue
            if i == 0:
                backtrack[i][j] = mark[2]
                s[i][j] = s[0][j-1] - gap
            elif j == 0:
                backtrack[i][j] = mark[1]
                s[i][j] = 0
            else:
                match = 1
                if a[i-1] != b[j-1]:
                    match = -1
                test = [s[i-1][j-1] + match,
                        s[i-1][j] - gap,
                        s[i][j-1] - gap]
                s[i][j] = max(test)
                
                index = test.index(max(test))
                backtrack[i][j] = mark[index]
    
    last_column = [x[-1] for x in s]
    max_score = max(last_column)
    n = last_column.index(max(last_column))
    
    align_a, align_b = '', ''
    while m != 0:
        if backtrack[n][m] == mark[0]:
            align_a += a[n-1]
            align_b += b[m-1]
            n -= 1
            m -= 1
        elif backtrack[n][m] == mark[1]:
            align_a += a[n-1]
            align_b += '-'
            n -= 1
        elif backtrack[n][m] == mark[2]:
            align_a += '-'
            align_b += b[m-1]
            m -= 1
        else:
            break
    
    return max_score, align_a[::-1], align_b[::-1]


def overlap_align(a, b, gap=2):
    n, m = len(a), len(b)
    s = [([0]*(m+1)) for i in range(n+1)]
    backtrack = [([' ']*(m+1)) for i in range(n+1)]
    mark = ['↘', '↓', '→']

    for i in range(n+1):
        for j in range(m+1):
            if i == 0 and j == 0:
                continue
            if i == 0:
                backtrack[i][j] = mark[2]
                s[i][j] = s[0][j-1] - gap
            elif j == 0:
                backtrack[i][j] = mark[1]
                s[i][j] = 0
            else:
                match = 1
                if a[i-1] != b[j-1]:
                    match = -gap
                test = [s[i-1][j-1] + match,
                        s[i-1][j] - gap,
                        s[i][j-1] - gap]
                s[i][j] = max(test)
                
                index = test.index(max(test))
                backtrack[i][j] = mark[index]
    
    last_row = s[-1]
    max_score = max(last_row)
    m = last_row.index(max(last_row))
    
    align_a, align_b = '', ''
    while m != 0:
        if backtrack[n][m] == mark[0]:
            align_a += a[n-1]
            align_b += b[m-1]
            n -= 1
            m -= 1
        elif backtrack[n][m] == mark[1]:
            align_a += a[n-1]
            align_b += '-'
            n -= 1
        elif backtrack[n][m] == mark[2]:
            align_a += '-'
            align_b += b[m-1]
            m -= 1
        else:
            break
    
    return max_score, align_a[::-1], align_b[::-1]


def alignment_score(str1, str2, match, mismatch, indel):
    score = []
    for a, b in zip(str1, str2):
        if a == b:
            score.append(match)
        elif a == '-' or b == '-':
            score.append(indel)
        else:
            score.append(mismatch)
    return sum(score)


def global_affine_gap_prot(a, b, gap_cost=11, gap_extension=1):
    with open('data\sm\BLOSUM62.pkl', 'rb') as input_file:
        score = pickle.load(input_file)
    
    n, m = len(a), len(b)
    s = [([0]*(m+1)) for i in range(n+1)]
    backtrack = [([' ']*(m+1)) for i in range(n+1)]
    mark = ['↘', '↓', '→']
    t, k = 1, 1
    
    for i in range(n+1):
        for j in range(m+1):
            # current affine penalty
            cost = (t*gap_cost) + ((k-1)*gap_extension)
            
            if i == 0 and j == 0:
                continue
            if i == 0:
                if j == 1:
                    s[i][j] = -gap_cost
                else:
                    backtrack[i][j] = mark[2]
                    s[i][j] = s[0][j-1] - gap_extension
            elif j == 0:
                if i == 1:
                    s[i][j] = -gap_cost
                else:
                    backtrack[i][j] = mark[1]
                    s[i][j] = s[i-1][0] - gap_extension
            else:
                test = [s[i-1][j-1] + score[a[i-1]][b[j-1]],
                        s[i-1][j] - cost,
                        s[i][j-1] - cost]
                s[i][j] = max(test)
                
                index = test.index(max(test))
                backtrack[i][j] = mark[index]
                
                # affine gap definition
                if index == 0:
                    t = 1
                    k = 1
                else:
                    t = 0
                    k += 1
                
    max_score = s[n][m]
    align_a, align_b = aligner(backtrack, a, b, n, m, mark)
    
    return max_score, align_a, align_b


def match_global_align(a, b):
    n, m = len(a), len(b)
    s = [([0]*(m+1)) for i in range(n+1)]
    backtrack = [([' ']*(m+1)) for i in range(n+1)]
    mark = ['↘', '↓', '→']

    for i in range(n+1):
        for j in range(m+1):
            if i == 0 and j == 0:
                continue
            if i == 0:
                backtrack[i][j] = mark[2]
                s[i][j] = s[0][j-1]
            elif j == 0:
                backtrack[i][j] = mark[1]
                s[i][j] = s[i-1][0]
            else:
                score = 0
                if a[i-1] == b[j-1]:
                    score = 1
                test = [s[i-1][j-1] + score,
                        s[i-1][j],
                        s[i][j-1]]
                s[i][j] = max(test)
                
                index = test.index(max(test))
                backtrack[i][j] = mark[index]
    
    max_score = s[n][m]
    
    align_a, align_b = aligner(backtrack, a, b, n, m, mark)
    
    return max_score, align_a, align_b


def mlcs(seq1, seq2, seq3):
    len1 = len(seq1)
    len2 = len(seq2)
    len3 = len(seq3)
    
    ### part 1: create backtrack matrix
    
    # backtrack matrix with blank initial values
    bt = [[['_'] * len3 for j in range(len2)]
             for i in range(len1)
         ]
    
    # score matrix with blank initial values
    s = [[[0] * (len3 + 1) for j in range(len2 + 1)]
            for i in range(len1 + 1)
        ]
    
    # fill in matrices cell by cell
    for i, ch1 in enumerate(seq1, start=1):
        for j, ch2 in enumerate(seq2, start=1):
            for k, ch3 in enumerate(seq3, start=1):
                match = 1 if (ch1 == ch2 == ch3) else 0
                cs = max(s[i][j-1][k-1],           # cs -current value
                         s[i-1][j][k-1],           # in score matrix
                         s[i-1][j-1][k],           #    =s[i][j][k]
                         s[i][j][k-1],
                         s[i][j-1][k],
                         s[i-1][j][k],
                         s[i-1][j-1][k-1] + match)
                
                if cs == s[i-1][j-1][k-1] + match:
                    cbt = 'd'                      # cbt -current value in
                elif cs == s[i][j-1][k-1]:         # backtrack matrix
                    cbt = '⇒'                     #    =bt[i-1][j-1][k-1]
                elif cs == s[i-1][j][k-1]:
                    cbt = '⇓'
                elif cs == s[i-1][j-1][k]:
                    cbt = '↘'
                elif cs == s[i][j][k-1]:
                    cbt = '⨂'
                elif cs == s[i][j-1][k]:
                    cbt = '→'
                elif cs == s[i-1][j][k]:
                    cbt = '↓'
                else:
                    raise ValueError('something wrong in calculation')
                    
                s[i][j][k] = cs
                bt[i-1][j-1][k-1] = cbt


    ### part 2: distribute gaps in the sequences and calculate score
    
    i, j, k = len1 - 1, len2 - 1, len3 - 1
    seq1_mod = seq2_mod = seq3_mod = ''   # modified versions of strings
    sc = 0                                # score
    
    while i != -1 and j != -1 and k != -1:
        ch = bt[i][j][k]
        if ch == '↓':
            seq1_mod += seq1[i]
            seq2_mod += '-'
            seq3_mod += '-'
            i -= 1
        elif ch == '→':
            seq1_mod += '-'
            seq2_mod += seq2[j]
            seq3_mod += '-'
            j -= 1
        elif ch == '⨂':
            seq1_mod += '-'
            seq2_mod += '-'
            seq3_mod += seq3[k]
            k -= 1
        elif ch == '↘':
            seq1_mod += seq1[i]
            seq2_mod += seq2[j]
            seq3_mod += '-'
            i -= 1
            j -= 1  
        elif ch == '⇓':
            seq1_mod += seq1[i]
            seq2_mod += '-'
            seq3_mod += seq3[k]
            i -= 1
            k -= 1
        elif ch == '⇒':
            seq1_mod += '-'
            seq2_mod += seq2[j]
            seq3_mod += seq3[k]
            j -= 1
            k -= 1
        elif ch == 'd':
            seq1_mod += seq1[i]
            seq2_mod += seq2[j]
            seq3_mod += seq3[k]
            if seq1[i] == seq2[j] == seq3[k]:
                sc += 1
            i -= 1
            j -= 1
            k -= 1
        else:
            raise ValueError(f'Unrecognize symbol in matrix: {ch}')

            
    # finalize lines after previous loop (simplest approach)
    while i >= 0 or j >= 0 or k >= 0:
        if (j==-1) and (k==-1):
            seq1_mod += seq1[i]
            seq2_mod += '-'
            seq3_mod += '-'
            i -= 1
        elif (i==-1) and (k==-1):
            seq1_mod += '-'
            seq2_mod += seq2[j]
            seq3_mod += '-'
            j -= 1
        elif (i==-1) and (j==-1):
            seq1_mod += '-'
            seq2_mod += '-'
            seq3_mod += seq3[k]
            k -= 1
        elif k==-1:
            seq1_mod += seq1[i]
            seq2_mod += seq2[j]
            seq3_mod += '-'
            i -= 1
            j -= 1  
        elif j==-1:
            seq1_mod += seq1[i]
            seq2_mod += '-'
            seq3_mod += seq3[k]
            i -= 1
            k -= 1
        elif i==-1:
            seq1_mod += '-'
            seq2_mod += seq2[j]
            seq3_mod += seq3[k]
            j -= 1
            k -= 1
        else:
            raise ValueError('Something wrong')

    return sc, seq1_mod[::-1], seq2_mod[::-1], seq3_mod[::-1]


def middle_node(in1, in2, top, bottom, left, right, score, gap):
    v = in1[top:bottom]
    w = in2[left:right]
    n, m = len(v), len(w)
    middle = math.floor(m/2)
    
    def to_middle(v, w):
        s = [None] * (n+1)
        for j in range(1, len(w)+1):
            if j == 1:
                p = list(range(0, +gap*(n+1), +gap))
            for i in range(n+1):
                if i == 0:
                    s[i] = j*gap
                else:
                    values = (p[i-1] + score[v[i-1]][w[j-1]],
                              p[i] + gap,
                              s[i-1] + gap)
                    
                    s[i] = max(values)
            p = copy.deepcopy(s)
        return s
    
    from_source = to_middle(v, w[:middle])
    to_sink = to_middle(v[::-1], w[middle:][::-1])
    sums = [x+y for x, y in zip(from_source, to_sink[::-1])]
    high = max(sums)
    middle_node = sums.index(high)
    
    return (middle_node, middle)


def greedy_sorting(P):
    R = list(range(1, len(P)+1))
    permutations = [P]
    
    def k_sort(P, k):
        if k in P:
            k_index = P.index(k)
        else:
            k_index = P.index(-k)
        flip = P[k-1:k_index+1][::-1]
        return P[:k-1] + [-i for i in flip] + P[k_index+1:]
    
    for pk, k in enumerate(R):
        if P[pk] == k:
            continue
        elif P[pk] == -k:
            P = P[:pk] + [k] + P[pk+1:]
            permutations.append(P)
        else:
            P = k_sort(P, k)
            permutations.append(P)
            if P[pk] == -k:
                P = P[:pk] + [k] + P[pk+1:]
                permutations.append(P)
        
    return permutations


def print_permutations(permutations):
    for P in permutations[1:]:
        print(' '.join(('+' if i > 0 else '') + str(i) for i in P))


def breakpoints(P):
    brkp, adj = 0, 0
    P = [0] + P + [len(P)+1]
    for i in range(len(P)-1):
        if P[i+1] - P[i] == 1:
            adj += 1
        else:
            brkp += 1
    return brkp, adj


def chromosome_to_cycle(chromosome):
    # print(chromosome)
    nodes = [None] * (2*(len(chromosome)))
    for j, i in enumerate(chromosome):
        if i > 0:
            nodes[(2*j)] = (2*i)-1
            nodes[(2*j)+1] = (2*i)
        else:
            nodes[(2*j)] = -(2*i)
            nodes[(2*j)+1] = -(2*i)-1
    return nodes


def cycle_to_chromosome(nodes):
    """
    Converts a list of directed edges in a chromosome into a list of synteny
    blocks and their orientation in that chromosome.

    Parameters
    ----------
    nodes : LIST
        List of nodes in a chromosome.

    Returns
    -------
    list
        DESCRIPTION.

    Example
    -------
    input:
        [1, 2, 4, 3, 6, 5]
        
    output:
        [1, -2, -3]

    """
    n = int(len(nodes)/2)
    chromosome = [None] * n
    for j in range(n):
        if nodes[2*j] < nodes[(2*j)+1]:
            chromosome[j] = nodes[(2*j) + 1]/2
        else:
            chromosome[j] = -nodes[(2*j)]/2
    return [int(i) for i in chromosome]


# def colored_edges(P):
#     edges = []
#     if isinstance(P[0], int):
#         temp = []
#         nodes = chromosome_to_cycle(P)
#         for j in range(len(P)):
#             temp.append((nodes[(2*j)-1], nodes[(2*j)]))
#         temp = temp[1:] + [temp[0]]
#         edges.append(temp)
#     else:
#         for chromosome in P:
#             temp = []
#             nodes = chromosome_to_cycle(chromosome)
#             for j in range(len(chromosome)):
#                 temp.append((nodes[(2*j)-1], nodes[(2*j)]))
#             temp = temp[1:] + [temp[0]]
#             edges.append(temp)
#     return edges


def colored_edges(P):
    edges = []
    for chromosome in P:
        temp = []
        nodes = chromosome_to_cycle(chromosome)
        for j in range(len(chromosome)):
            temp.append((nodes[(2*j)-1], nodes[(2*j)]))
        temp = temp[1:] + [temp[0]]
        edges.append(temp)
    return [item for ele in edges for item in ele]


def graph_to_genome(genome_graph):
    """
    Converts a list of ordered undirected nodes in a genome graph to a genome
    of synteny blocks separated into chromosomes (if applicable).
    
    Dependencies
    ------------
        cycle_finder(graph)
        cycle_to_genome(cycle)

    Parameters
    ----------
    genome_graph : LIST of tuples
        An ordered list of undirected nodes in genome graph in cyclical format.

    Returns
    -------
    P : LIST of lists
        List of lists where each list represents a chromosome and each
        chromosome represents an ordering of synteny blocks within that
        chromosome.
        
    Example
    -------
    input:
        [(2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8)]
        
    output:
        [[1, -2, -3], [-4, 5, -6]]
        
    """
    P = []
    cycles = cycle_finder(genome_graph)
    # for each cycle in the graph
    for i in cycles:
        # convert cycle to directed edges of synteny blocks
        cycle = list(itertools.chain(*cycles[i]))
        chromosome = cycle_to_chromosome([cycle[-1]] + cycle[:-1])
        P.append(chromosome)
    return P


def cycle_finder(graph):
    """
    Splits an ordered genome graph into its independent cycles by taking the
    difference between the outgoing nodes value and the incoming value of the
    next nodes. if this value-the difference- is greater than 1, a new cycle
    is found.

    Parameters
    ----------
    graph : LIST
        List of ordered nodes in genome_graph.

    Returns
    -------
    cycles : DICT
        Dictionary where cycles[i].values() contains a list of nodes (tuples)
        in cycle i. NB, cycles are 1-indexed, so the first cycle of graph
        will be cycles[1] NOT cycles[0].
        
    Example
    -------
    input:
        [(2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8)]
    
    output:
        {1: [(2, 4), (3, 6), (5, 1)], 2: [(7, 9), (10, 12), (11, 8)]}

    """
    nodes = copy.deepcopy(graph)
    cycles = {}
    k = 0
    while nodes:
        k += 1
        node = nodes[0]
        cycle = []
        while node not in cycle:
            cycle.append(node)
            if node[1]%2 == 0:
                next_node = [i for i in nodes if i[0] == node[1]-1]
            else:
                next_node = [i for i in nodes if i[0] == node[1]+1]
            nodes.pop(nodes.index(node))
            if next_node:
                node = next_node[0]
            else:
                cycles[k] = cycle
                break
    return cycles


def two_break(graph, i1, i2, i3, i4):
    """
    Rearranges nodes in graph based on a two break, changing nodes (i1, i2)
    and (i3, i4) to (i1, i3) and (i2, i4), respectively. Retaining proper
    oredering and orientation throughout

    Parameters
    ----------
    graph : LIST of tuples
        List of tuples representing undirected nodes in the genome graph.
    i1 : INT
    i2 : INT
    i3 : INT
    i4 : INT
        Where (i1, i2) and (i3, i4) are in the genome graph

    Returns
    -------
    graph : LIST of nodes
        New genome graph after 2-break has been implemented at (i1, i2) and 
        (i3, i4).
        
    Example
    -------
    input:
        [(2, 4), (3, 8), (7, 5), (6, 1)], 1, 6, 3, 8
        
    output:
        [(2, 4), (3, 8), (7, 5), (6, 8)]

    """
    graph = copy.deepcopy(graph)
    if (i1, i2) in graph:
        index = graph.index((i1, i2))
        graph[index] = (i4, i2)
    elif (i2, i1) in graph:
        index = graph.index((i2, i1))
        graph[index] = (i2, i4)
    if (i3, i4) in graph:
        index = graph.index((i3, i4))
        graph[index] = (i3, i1)
    elif (i4, i3) in graph:
        index = graph.index((i4, i3))
        graph[index] = (i1, i3)
    return graph


def genome_two_break(P, i1, i2, i3, i4):
    """
    Function applying a two break on genome P, changing nodes (i1, i2)
    and (i3, i4) to (i1, i3) and (i2, i4), respectively.

    Parameters
    ----------
    P : LIST (of lists-if multiple chromosomes)
        Genome P made of orientated synteny blocks on potentially multiple
        chromosomes.
    i1 : INT
    i2 : INT
    i3 : INT
    i4 : INT
        Where (i1, i2) and (i3, i4) are in the genome graph

    Returns
    -------
    P : LIST (of lists)
        New genome P after 2-break.

    Example
    -------
    input:
        [1, -2, -4, 3], 1, 6, 3, 8
        
    output:
        [[1, -2], [-4, 3]]
        
    """
    graph = colored_edges(P)
    graph = two_break(graph, i1, i2, i3, i4)
    P = graph_to_genome(graph)
    return P


def two_break_distance(P, Q):
    red = colored_edges(P)
    block = len(red)
    blue = colored_edges(Q)
    G = red + blue
    cycles = len(cycle_finder(G))
    return block - cycles


def common(a,b): 
    c = [value for value in a if value in b] 
    return c


def shortest_rearrangement(P, Q):
    out = P
    red_edges = colored_edges(P)
    blue_edges = colored_edges(Q)
    breakpoint_graph = red_edges + blue_edges
    cycles = cycle_finder(breakpoint_graph)
    while len(cycles) < len(red_edges): # i.e. non-trivial cycles
        non_trivial_cycles = [cycles[i] for i in cycles if len(cycles[i]) > 2]
        non_trivial_blue_edges = common(blue_edges, *non_trivial_cycles)
        arbitrary_edges = (non_trivial_blue_edges)
        for edge in arbitrary_edges:
            red_edge_1 = [i for i in red_edges if i[1] == edge[0]]
            red_edge_2 = [i for i in red_edges if i[0] == edge[1]]
            if red_edge_1 and red_edge_2:
                break
        new_edge_1 = edge
        new_edge_2 = (red_edge_2[0][1], red_edge_1[0][0])
        red_edges[red_edges.index(*red_edge_1)] = new_edge_1
        red_edges[red_edges.index(*red_edge_2)] = new_edge_2
        breakpoint_graph = red_edges + blue_edges
        i1 = red_edge_1[0][0]
        i2 = red_edge_1[0][1]
        i3 = red_edge_2[0][0]
        i4 = red_edge_2[0][1]
        out.append(P)
        P = (genome_two_break(P, i1, i2, i3, i4))
    return out


def shared_kmers(k, a, b):
    out = []
    kmers = I.all_kmers(a, k)
    rcs = I.all_kmers(I.reverse_complement(a), k)[::-1]
    for i, (kmer, rc) in enumerate(zip(kmers, rcs)):
        ys = []
        if kmer in b:
            ys.append([n for n in range(len(b)) if b.find(kmer, n) == n])
        if rc in b:
            ys.append([n for n in range(len(b)) if b.find(rc, n) == n])
        if ys:
            for line in ys:
                out.append([(i, y) for y in line])
            del(ys)
        else:
            continue
    return [item for ele in out for item in ele]



# =============================================================================
# Testing
# =============================================================================
if __name__ == "__main__":
    my.flash()
    
    a = 'TGCCCCGGTGGTGAG'
    b = 'AAGGTCGCACCTCGT'
    print(len(shared_kmers(3, a, b)))
    # my.list_writer(out)
    
    # raw = my.reader('input.txt')
    # genome = []
    
    # for line in raw:
    #     line = line[1:-1]
    #     line = line.split(')(')
    #     gen = []
    #     for chromosome in line:
    #         temp = chromosome.split()
    #         out = [int(i) for i in temp]
    #         gen.append(out)
            
    #     genome.append(gen)
    
    # P = genome[0]
    # Q = genome[1]
    
    # # out = (two_break_distance(P, Q))
    # # out = (colored_edges(Q))
    
    # out = (shortest_rearrangement(P, Q))
    
    # print(out)
    
    # out = (genome_two_break(edges, 1, 6, 3, 8))
    # out = (genome_two_break(edges, 24, 33, 125, 4))
    
    # submit = ''
    # for line in out:
    #     submit += '('
    #     for i, ele in enumerate(line):
    #         if ele < 0:
    #             submit += str(int(ele))
    #         else:
    #             submit += '+' + str(ele)
    #         if i == len(line)-1:
    #             submit += ')'
    #         else:
    #             submit += ' '
    
    # print(submit)
    
    # cycle = [(2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8)]
    # print(graph_to_genome(cycle))
    
    # data = my.reader('input.txt')
    # data = data[0][1:-1]
    # data = data.split()
    # data = [int(i) for i in data]
    
    # out = (chromosome_to_cycle(data))
    
    # print(out)
    
    
    
    