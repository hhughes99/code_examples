import helper as my
import networkx as nx

# =============================================================================
# SECTION 1
# =============================================================================
def composition(string, k):
    n = len(string)
    kmers = []
    for i in range(n-(k-1)):
        kmers.append(string[i: i+k])
    kmers.sort()
    return kmers

def path_to_genome(path):
    n = len(path)-1
    genome = ''
    for i, string in enumerate(path):
        if i == n:
            genome += string
        else:
            genome += string[0]
    return genome

def overlap(kmers):
    out = {}
    for kmer in kmers:
        for comp in kmers:
            if kmer[1:] == comp[:-1]:
                if kmer in out:
                    out[kmer].append(comp)
                else:
                    out[kmer] = comp.split()
    return out

def all_kmers(string, k):
    kmers = []
    for i in range(len(string)-(k-1)):
        kmers.append(string[i: i + k])
    return kmers

def debruijn(k, text):
    out = {}
    edges = all_kmers(text, k)
    nodes = [x[:k-1] for x in edges]
    nodes.append(text[-k+1:])
    for node in set(nodes):
        for edge in edges:
            prefix = edge[:-1]
            if node == prefix:
                if node in out:
                    out[node].append(edge[1:])
                else:
                    out[node] = edge[1:].split()
    return out

def debruijn_from_kmers(kmers):
    d = {}
    for kmer in kmers:
        if kmer[:-1] not in d:
            d[kmer[:-1]] = [kmer[1:]]
        else:
            d[kmer[:-1]].append(kmer[1:])
    return d

def dict_to_elist(d):
    out = []
    for entry in d:
        for value in d[entry]:
            out.append((entry, value))
    return out

def adj_reader(strings):
    d = {}
    for string in strings:
        string = string.split(' -> ')
        d[string[0]] = string[1].split(',')
    return d

def eulerian_cycle(d):
    elist = dict_to_elist(d)
    G = nx.DiGraph()
    G.add_edges_from(elist)
    cycle = [u for u, v in nx.eulerian_circuit(G)]
    cycle.append(cycle[0])
    return cycle

def eulerian_path(d):
    elist = dict_to_elist(d)
    G = nx.DiGraph()
    G.add_edges_from(elist)
    # G = nx.eulerize(G)
    path = list(nx.eulerian_path(G))
    out = [u for u, v in path]
    out.append(path[-1][-1])
    return out

def cycle_format(cycle):
    return '->'.join(cycle)

def string_reconstruction(patterns):
    db = debruijn_from_kmers(patterns)
    path = eulerian_path(db)
    out = path_to_genome(path)
    return out

def cycle_reconstruction(patterns):
    k = len(patterns[0])
    db = debruijn_from_kmers(patterns)
    path = eulerian_cycle(db)
    out = path_to_genome(path)
    return out[:-k+1]

def hamming(seq1, seq2):
    return sum([1 for x, y in zip(seq1, seq2) if x != y])

def binary_neighbors(pattern, d):
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return {'1', '0'}
    neighborhood = set()
    suffix_neighbors = binary_neighbors(pattern[1:], d)
    for text in suffix_neighbors:
        if hamming(pattern[1:], text) < d:
            for x in '01':
                neighborhood.add(x + text)
        else:
            neighborhood.add(pattern[0] + text)
    return neighborhood

def binary_strings(k):
    pattern = '0' * k
    return list(binary_neighbors(pattern, k))

def paired_composition(string, k, d):
    n = len(string)
    pairs = []
    for i in range(n-(2*k)-d+1):
        pairs.append((string[i:i+k], string[i+k+d:i+(2*k)+d]))
    pairs.sort()
    return pairs

def genome_from_gapped_patterns(gapped_patterns, k, d):
    first = [x[0] for x in gapped_patterns]
    second = [x[1] for x in gapped_patterns]
    prefix = path_to_genome(first)
    suffix = path_to_genome(second)
    for i in range((k+d), len(prefix)):
        if prefix[i] != suffix[i-k-d]:
            return 'No string spelled by gapped patterns'
    return prefix[:k+d] + suffix

def debruijn_from_pairs(pairs):
    d = {}
    for pair in pairs:
        prefix = pair[0][:-1] + '|' + pair[1][:-1]
        suffix = pair[0][1:] + '|' + pair[1][1:]
        if prefix not in d:
            d[prefix] = suffix
        else:
            d[prefix] += (suffix)
    return d
        


# =============================================================================
# Testing
# =============================================================================
if __name__ == "__main__":
    my.flash()
    
    # raw = my.reader('dataset_203_7.txt')
    # pairs = 'AG|AG → GC|GC → CA|CT → AG|TG → GC|GC → CT|CT → TG|TG → GC|GC → CT|CA'
    
    pairs = my.reader('input.txt')
    
    data = []
    for pair in pairs:
        data.append(pair.split('|'))
    
    d = (debruijn_from_pairs(data))
        
    eulerian_cycle(d)
    
    
    
    
    
    # # eulerian cycle/path output pipeline
    # raw = my.reader('dataset_203_6.txt')
    # d = adj_reader(raw)
    # path = (eulerian_path(d))
    # out = (cycle_format(path))
    # my.writer(out)