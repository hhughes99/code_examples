import helper as my
import networkx as nx
import numpy as np
import copy

# =============================================================================
# SECTION 1
# =============================================================================
def distance_matrix(G):
    leaves = [x for x in G.nodes() if G.out_degree(x)==1 and G.in_degree(x)==1]
    n = len(leaves)
    D = [([None]*(n)) for i in range(n)]
    for i, source in enumerate(leaves):
        for j, target in enumerate(leaves):
            if i == j:
                D[i][j] = 0
            else:
                D[i][j] = nx.shortest_path_length(G, source, target, weight='weight')
    return D


def limb_length(j, D):
    '''
    Given an additive matrix D and a leaf j, LimbLength(j) is equal to the 
    minimum value of:
    
    (D_{i,j} + D_{j,k} - D_{i,k})/2
    
    over all leaves i and k.

    Parameters
    ----------
    j : STRING
        Leaf in tree(D).
    D : LIST of LISTS
        Distance matrix.

    Returns
    -------
    None.

    '''
    n = len(D)
    stack = []
    for i in range(n):
        for k in range(n):
            if j == i or j == k:
                continue
            else:
                stack.append((D[i][j] + D[j][k] - D[i][k])/2)
    return int(min(stack))


def additive_phylogeny(D):
    n = len(D)-1
    if n == 2:
        return T.add_edge(0, 1, weight = D[0][1])
    
    limb = limb_length(n, D)
    
    for j in range(n):
        D[j][n] -= limb
        D[n][j] = D[j][n]
    
    for i in range(n + 1):
        for k in range(n + 1):
            if i == k:
                continue
            elif D[i][k] == D[i][n] + D[n][k]:
                break
    
    x = D[i][n]
    D = [row[:-1] for row in D[:-1]]
    
    T = additive_phylogeny(D)
    
    nx.draw(T)
    
    print(nx.dijkstra_path(T, i, k))
        
    return T
    

def discrepancy(T, D):
    stack = []
    for i in range(len(T)):
        for j in range(len(T[0])):
            stack.append((T[i][j]-D[i][j])**2)
    return sum(stack)/2


def lowest_cell(table):
    min_cell = float("inf")
    x, y = -1, -1
    for i in range(len(table)):
        for j in range(len(table[i])):
            if i == j:
                continue
            elif table[i][j] < min_cell:
                min_cell = table[i][j]
                x, y = i, j
    return x, y


class edge:
    def __init__(self, start, end, weight):
        self.start = start
        self.end = end
        self.weight = weight


def path_length(T, i):
    G = copy.deepcopy(T)
    G = [edge for edge in G if edge.end < edge.start]
    if i not in [i.start for i in G]:
        return 0
    
    d, i_stack = 0, []
    while True:
        starts = [i.start for i in G]
        if i not in starts:
            break
        else:
            find = starts.index(i)
            d += G[find].weight
            i = G[find].end
            i_stack.append(i)
            G.pop(find)
        
    return d


def upgma(D):
    T = []
    n = len(D)
    labels = list(range(n))
    clusters = {i:str(i) for i in labels}
    
    while len(D) > 1:
        
        j, i = lowest_cell(D)
        clusters[n] = clusters[i] + clusters[j]
        
        d = D[i][j]
        
        ic = path_length(T, labels[i])
        jc = path_length(T, labels[j])
        
        T.append(edge(n, labels[i], (d/2)-ic))
        T.append(edge(n, labels[j], (d/2)-jc))
        
        children = [i.end for i in T if i.start == n]
        
        if labels[i] in children:
            T.append(edge(labels[i], n, (d/2)-ic))
        if labels[j] in children:
            T.append(edge(labels[j], n, (d/2)-jc))
        
        new_distances = []
        for k in range(len(D)):
            if k == i or k == j:
                continue
            else:
                cluster_i = clusters[labels[i]]
                cluster_j = clusters[labels[j]]
                num = D[k][i]*len(cluster_i) + D[k][j]*len(cluster_j)
                den = len(cluster_i) + len(cluster_j)
                new_distances.append(num/den)   
        new_distances.append(0)
        
        labels.pop(i); labels.pop(j)
        labels.append(n)
        
        D = [[ele for index, ele in enumerate(r) if index != i and index != j] 
                 for index, r in enumerate(D) if index != i and index != j]
        
        for k, row in enumerate(D):
            row.append(new_distances[k])
        D.append(new_distances)
        
        n += 1
    
    T.sort(key=lambda x: (x.start, x.end))
    
    return T


def d_star(D):
    D2 = copy.deepcopy(D)
    n = len(D)
    
    d = [sum(i) for i in D]
    
    for i in range(n):
        for j in range(n):
            if j == i:
                D2[i][j] = 0
            else:
                D2[i][j] = (n-2)*D[i][j] - d[i] - d[j]
    
    return D2, d


def nj(D):
    T = []
    n = len(D)
    labels = list(range(n))
    clusters = {i:str(i) for i in labels}
    
    while len(D) > 1:
        
        j, i = lowest_cell(d_star(D)[0])
        clusters[n] = clusters[i] + clusters[j]
        
        d = D[i][j]
        
        ic = path_length(T, labels[i])
        jc = path_length(T, labels[j])
        
        T.append(edge(n, labels[i], (d/2)-ic))
        T.append(edge(n, labels[j], (d/2)-jc))
        
        children = [i.end for i in T if i.start == n]
        
        if labels[i] in children:
            T.append(edge(labels[i], n, (d/2)-ic))
        if labels[j] in children:
            T.append(edge(labels[j], n, (d/2)-jc))
        
        new_distances = []
        for k in range(len(D)):
            if k == i or k == j:
                continue
            else:
                cluster_i = clusters[labels[i]]
                cluster_j = clusters[labels[j]]
                num = D[k][i]*len(cluster_i) + D[k][j]*len(cluster_j)
                den = len(cluster_i) + len(cluster_j)
                new_distances.append(num/den)   
        new_distances.append(0)
        
        labels.pop(i); labels.pop(j)
        labels.append(n)
        
        D = [[ele for index, ele in enumerate(r) if index != i and index != j] 
                 for index, r in enumerate(D) if index != i and index != j]
        
        for k, row in enumerate(D):
            row.append(new_distances[k])
        D.append(new_distances)
        
        n += 1
    
    T.sort(key=lambda x: (x.start, x.end))
    
    return T




# =============================================================================
# Testing
# =============================================================================
if __name__ == "__main__":
    my.flash()
    
    raw = my.reader()
    data = []
    for line in raw:
        data.append([int(i) for i in line.split(' ')])
    del(line); del(raw)
    
    # out = (nj(data))
    
    # out = (upgma(data))
    
    # print('\n')
    # for line in out:
    #     print(str(line.start) + '->' + str(line.end) 
    #           + ':' + '%0.03f' % line.weight)
    
    T = [[0, 14, 17, 11],
         [14, 0, 21, 15],
         [17, 21, 0, 18],
         [11, 15, 18, 0]]
    
    print(discrepancy(T, data))
    
    D = [[0, 20, 9, 11],
         [20, 0, 17, 11],
         [9, 17, 0, 8],
         [11, 11, 8, 0]]
    
    print(d_star(D)[0][1][0])
    
    print(limb_length(1, D))
    
    
    
    
    