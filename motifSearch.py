import helper as my
import collections
import random
import math
from copy import copy

# =============================================================================
# SECTION 1
# =============================================================================
def uraciler(a):
    u = ''
    for letter in a:
        if letter == 'T':
            u += 'U'
        else:
            u += letter
    return u

def pattern_count(text, pattern):
    count = 0
    for i in range(len(text)-(len(pattern)-1)): # use of -1? # denoted REF
        if text[i: i+len(pattern)] == pattern:
            count += 1
    return count

def frequent_words(text, k):
    frequent_patterns = set()
    count = []
    for i in range(len(text)-k): #REF
        pattern = text[i: i+k]
        count.append(pattern_count(text, pattern))
    max_count = max(count)
    for i in range(len(text)-k): #REF
        if count[i] == max_count:
            frequent_patterns.add(text[i: i+k])
    return frequent_patterns       

def frequency_table(text, k):
    freq_map = {}
    n = len(text)
    for i in range(n-k): # REF
        pattern = text[i: i+k]
        if pattern not in freq_map:
            freq_map[pattern] = 1
        else:
            freq_map[pattern] += 1
    return freq_map

def better_frequent_words(text, k):
    frequent_patterns = []
    freq_map = frequency_table(text, k)
    max_freq = max(freq_map.values())
    for pattern in freq_map.keys():
        if freq_map[pattern] == max_freq:
            frequent_patterns.append(pattern)
    return frequent_patterns

def reverse_complement(sequence):
    complement = ''
    for base in sequence.upper():
        if base == 'A':
            complement += 'T'
        elif base == 'T':
            complement += 'A'
        elif base == 'C':
            complement += 'G'
        else:
            complement += 'C'
    return complement[::-1]

def pattern_starts(text, pattern):
    starts = []
    for i in range(len(text)-(len(pattern)-1)): # REF
        if text[i: (i+len(pattern))] == str(pattern):
            starts.append(i)
    return starts

def find_clumps(text, k, L, t):
    pattern = set()
    n = len(text)
    for i in range(n-L): # REF
        window = text[i: i+L]
        freq_map = frequency_table(window, k)
        for key in freq_map.keys():
            if freq_map[key] >= t:
                pattern.add(key)
    return pattern

# =============================================================================
# SECTION 2
# =============================================================================
def skew(text):
    score = 0
    skew = [0]
    for base in text.lower():
        if base == 'c':
            score -= 1
        elif base == 'g':
            score += 1
        skew.append(score)
    return skew

def find_index(lst, value):
    index = []
    for x, i in enumerate(lst):
        if i == value:
            index.append(x)
    return value, index

def ori_finder(text):
    s = skew(text)
    return find_index(s, min(s))[1]

def hamming_distance(seq1, seq2):
    diff = 0
    for a, b in zip(seq1, seq2):
        if a != b:
            diff += 1
    return diff

def approx_pattern_starts(text, pattern, d):
    starts =[]
    for i in range(len(text)-(len(pattern)-1)): #REF
        if hamming_distance(text[i: (i+len(pattern))], pattern) <= d:
            starts.append(i)
    return starts

def approx_pattern_count(text, pattern, d):
    count = 0
    for i in range(len(text)-(len(pattern)-1)): # REF
        if hamming_distance(text[i: (i+len(pattern))], pattern) <= d:
            count += 1
    return count

def immediate_neighbors(pattern):
    neighborhood = []
    for i in range(len(pattern)):
        symbol = pattern[i]
        for x in 'ATCG':
            if symbol != x:
                neighbor = pattern[:i] + x + pattern[i+1:]
                neighborhood.append(neighbor)
    return neighborhood

def neighbors(pattern, d):
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return {'A', 'C', 'G', 'T'}
    neighborhood = set()
    suffix_neighbors = neighbors(pattern[1:], d)
    for text in suffix_neighbors:
        if hamming_distance(pattern[1:], text) < d:
            for x in 'ATCG':
                neighborhood.add(x + text)
        else:
            neighborhood.add(pattern[0] + text)
    return neighborhood

def all_strings(k):
    pattern = 'A' * k
    return neighbors(pattern, k)
    
def frequent_mismatches(text, k, d):
    patterns = set()
    freq_map = {}
    n = len(text)
    for i in range(n-k): # REF
        pattern = text[i: i+k]
        neighborhood = neighbors(pattern, d)
        for neighbor in neighborhood:
            if neighbor not in freq_map:
                freq_map[neighbor] = 1
            else:
                freq_map[neighbor] += 1
    m = max(freq_map.values())
    for pattern in freq_map.keys():
        if freq_map[pattern] == m:
            patterns.add(pattern)
    return patterns
    
def frequent_mismatches_with_rc(text, k, d):
    patterns = set()
    freq_map = {}
    n = len(text)
    for i in range(n-k): # REF
        for pattern in [text[i: i+k], reverse_complement(text[i: i+k])]:
            neighborhood = neighbors(pattern, d)
            for neighbor in neighborhood:
                if neighbor not in freq_map:
                    freq_map[neighbor] = 1
                elif reverse_complement(neighbor) not in freq_map:
                    freq_map[neighbor] = 1
                else:
                    freq_map[neighbor] += 1
    m = max(freq_map.values())
    for pattern in freq_map.keys():
        if freq_map[pattern] == m:
            patterns.add(pattern)
    return patterns

# =============================================================================
# SECTION 3
# =============================================================================
def motif_enumeration(dna, k, d):
    kmer_list = [set() for _ in dna]
    for i, seq in enumerate(dna):
        for j in range(len(seq)-(k-1)):
            k_mer = seq[j: j+k]
            
            kmer_list[i].update(neighbors(k_mer, d))
    return set.intersection(*kmer_list)

def count_motifs(motifs):
    profile = {'A': [], 'C': [], 'G': [], 'T': []}
    transpose = [''.join(s) for s in zip(*motifs)]
    for row in transpose:
        for base in 'ACGT':
            profile[base].append(row.count(base))
    return profile

def profile_motifs(motifs):
    profile = {'A': [], 'C': [], 'G': [], 'T': []}
    if not isinstance(motifs, list):
        motifs = list(motifs.split(" "))
    transpose = [''.join(s) for s in zip(*motifs)]
    for row in transpose:
        for base in 'ACGT':
            profile[base].append((row.count(base))/len(row))
    return profile

def consensus(motifs):
    consensus = ''
    transpose = [''.join(s) for s in zip(*motifs)]
    for row in transpose:
        base = collections.Counter(row).most_common(1)[0][0]
        consensus += base
    return consensus

def calc_entropy(prob):
    return round(-sum([x*math.log2(x) for x in prob if x > 0]), 4)

def entropy(motifs):
    entropy = []
    transpose = [''.join(s) for s in zip(*motifs)]
    for row in transpose:
        lst = []
        for base in 'acgt':
            lst.append(row.lower().count(base)/len(row))
        entropy.append(calc_entropy(lst))
    return sum(entropy), entropy

# all_kmers = lambda k: neighbors('a' * k, k)

def d_pattern_motifs(pattern, motifs):
    d = 0
    for motif in motifs:
        d += hamming_distance(pattern, motif)
    return d

def d_pattern_strings(pattern, dna):
    k = len(pattern)
    d = 0
    for string in dna:
        hamming_d = float('inf')
        for i in range(len(string)-(k-1)):
            kmer = string[i: i+k]
            if hamming_d > hamming_distance(kmer, pattern):
                hamming_d = hamming_distance(kmer, pattern)
        d += hamming_d
    return d

def median_string(dna, k):
    d = float('inf')
    kmers = all_strings(k)
    median = ''
    for kmer in kmers:
        if d > d_pattern_strings(kmer, dna):
            d = d_pattern_strings(kmer, dna)
            median = kmer
    return median

def all_kmers(string, k):
    kmers = []
    for i in range(len(string)-(k-1)):
        kmers.append(string[i: i + k])
    return kmers

def greedy_motif_search(dna, k):
    t = len(dna)
    best_motifs = [x[0:k] for x in dna]
    for kmer in all_kmers(dna[0], k):
        motifs = [kmer]
        for i in range(1, t):
            profile = profile_motifs(motifs)
            motifs.append(profile_most_probable(dna[i], profile, k))
        if score_motifs(motifs) < score_motifs(best_motifs):
            best_motifs = motifs
    return best_motifs

def greedy_motif_search_ii(dna, k):
    t = len(dna)
    best_motifs = [x[0:k] for x in dna]
    for kmer in all_kmers(dna[0], k):
        motifs = [kmer]
        for i in range(1, t):
            profile = laplace_profile(motifs)
            motifs.append(profile_most_probable(dna[i], profile, k))
        if score_motifs(motifs) < score_motifs(best_motifs):
            best_motifs = motifs
    return best_motifs

# =============================================================================
# SECTION 4
# =============================================================================
def random_kmers(dna, k):
    kmers = []
    for string in dna:
        i = random.randint(0, (len(string)-(k)))
        kmers.append(string[i: i+k])
    return kmers

def score_motifs(motifs):
    score = []
    transpose = [''.join(s) for s in zip(*motifs)]
    for row in transpose:
        freq = collections.Counter(row).most_common(1)[0][1]
        score.append(len(row)-freq)
    return sum(score)

def laplace_profile(motifs):
    profile = {'A': [], 'C': [], 'G': [], 'T': []}
    transpose = [''.join(s) for s in zip(*motifs)]
    for row in transpose:
        for base in 'ACGT':
            value = (row.count(base)+1)/(len(row)+4)
            profile[base].append(value)
    return profile

def pr(motif, profile):
    lst = []
    for i, base in enumerate(motif):
        if not i:
            lst = profile[base][i]
        else:
            lst *= profile[base][i]
    return lst

def profile_most_probable(string, profile, k):
    p = 0
    kmer = string[:k]
    for i in range(len(string)-(k-1)):
        motif = string[i: i+k]
        if p < pr(motif, profile):
            p = pr(motif, profile)
            kmer = motif
    return kmer

def profile_random_string(string, profile, k):
    p = []
    motif = []
    for i in range(len(string)-(k-1)):
        motif.append(string[i: i+k])
        p.append(pr(motif[i], profile))
    return random.choices(motif, p)[0]
    
def profile_most_probable_strings(strings, profile, k):
    motifs = []
    for string in strings:
        motifs.append(profile_most_probable(string, profile, k))
    return motifs

def random_motif_search(dna, k):
    motifs = random_kmers(dna, k)
    init_score = score_motifs(motifs)
    while True:
        profile = laplace_profile(motifs)
        motifs = profile_most_probable_strings(dna, profile, k)
        score = score_motifs(motifs)
        if score < init_score:
            init_score = score
        else:
            return init_score, motifs

def monte_carlo_motif_search(dna, k, t):
    for i in range(t):
        iteration = random_motif_search(dna, k)
        if not i:
            out = iteration
        if iteration[0] < out[0]:
            out = iteration
    return out

def gibbs_random(lst, p):
    C = sum(p)
    p = [x/C for x in p]
    out = random.choices(lst, p)
    return out[0]

def gibbs_sampler(dna, k, t, N):
    motifs = random_kmers(dna, k)
    best_motifs = copy(motifs)
    for j in range(N):
        i = random.randint(0, t-1)
        profile = laplace_profile(motifs[:i] + motifs[i+1:])
        motifs[i] = profile_random_string(dna[i], profile, k)
        if score_motifs(motifs) < score_motifs(best_motifs):
            best_motifs = copy(motifs)
    return score_motifs(best_motifs), best_motifs

def gibbs(dna, k, t, N, starts):
    for i in range(starts):
        iteration = gibbs_sampler(dna, k, t, N)
        if not i:
            out = iteration
        if iteration[0] < out[0]:
            out = iteration
    return out

# =============================================================================
# Testing
# =============================================================================
if __name__ == "__main__":
    my.flash()
    
    
    
    
    
    
    
    
    