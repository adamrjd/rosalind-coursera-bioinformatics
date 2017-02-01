from Parsing import MassTable, codons
from os import getcwd

lookup = MassTable()
rbase = 'ACGU'
base = 'ACGT'

global lookup
global rbase
global base

# ... SCRIPTS ......................................................
# Manipulating graphs


def circularize(element):
    if len(element) == 1:
        return element
    else:
        return [element[i:] + element[:i] for i in range(0, len(element))]


def decircularize(peptides):
    peptides = sorted(peptides)
    for i in range(0, len(peptides) / len(peptides[0])):
        cyclopeptide = peptides[i]
        removable = circularize(cyclopeptide)
        [peptides.remove(pep) for pep in removable
         if pep != cyclopeptide and pep in peptides]
    return peptides


def assemble(paths):
    reads = []
    for path in paths:
        read = path[0] + "".join([kmer[len(kmer) - 1:] for kmer in path[1:]])
        reads.append(read)
    return reads

# DNA manipulation numerically


def skew(dna):
    s = [0]
    skew_dict = {'G': 1, 'C': -1, 'A': 0, 'T': 0}
    for i, let in enumerate(dna):
        s.append(s[i] + skew_dict[let])
    return s


def hammingdistance(p, q):
    if len(p) != len(q):
        return 0
    else:
        return sum([1 for i in range(len(p)) if p[i] != q[i]])


def numbertopattern(num, k):
    return "".join(
        [numnuc[int(num / 4**(k - i)) % 4] for i in range(1, k + 1)])


def patterntonumber(path):
    l = len(path)
    return sum([nucnum[path[i]] * 4**(l - i - 1) for i in range(0, l)])

# RNA scripts


def reversecomplement(dna):
    return "".join([pairs[letter] for i, letter in enumerate(dna[::-1])])


def reversecomplementrna(rna):
    return "".join([rpairs[letter] for i, letter in enumerate(rna[::-1])])


def translate(rna):
    return "".join([codons[rna[i:i + 3]] for i in range(0, len(rna) - 2, 3)
                    if codons[rna[i:i + 3]] != '#'])


def transcribe(dna):
    return "".join([dnarna[let] for let in dna])


def reversetranscribe(rna):
    return "".join([dnarna[let] for let in rna])


def shared_kmers(k, dna1, dna2):
    kmers_dict = dict()
    shared_kmers_lst = list()
    enumerate_dna = lambda dna: enumerate(
        [dna[_:_ + k] for _ in range(len(dna) - k + 1)])
    for i, kmer in enumerate_dna(dna1):
        if kmer in kmers_dict:
            kmers_dict[kmer].append(i)
        else:
            kmers_dict[kmer] = list([i])
    for j, kmer in enumerate_dna(dna2):
        if kmer in kmers_dict:
            for x in kmers_dict[kmer]:
                shared_kmers_lst.append((x, j))
        elif reversecomplement(kmer) in kmers_dict:
            for x in kmers_dict[reversecomplement(kmer)]:
                shared_kmers_lst.append((x, j))
    return list(sorted(shared_kmers_lst, key=lambda l: l[1]))

# ... ROSALIND SOLUTIONS ....................................


def neighbors(Pattern, d):
    k = len(Pattern)
    letters = ["A", "C", "G", "T"]
    if d == 0:
        return [Pattern]

    if k == 1:
        return letters

    neighborhood = []
    suffixneighbors = neighbors(Pattern[1:len(Pattern) + 1], d)
    for i in range(0, len(suffixneighbors)):
        if hammingdistance(Pattern[1:len(Pattern) + 1],
                           suffixneighbors[i]) < d:
            for j in range(0, 4):
                neighborhood.append(letters[j] + suffixneighbors[i])
        else:
            neighborhood.append(Pattern[0:1] + suffixneighbors[i])
    return neighborhood


def approximate_pattern_count(Text, pattern, d):
    count = 0
    for i in range(0, len(Text) - len(pattern) + 1):
        pat = Text[i:i + len(pattern)]
        if hammingdistance(pat, pattern) <= d:
            count += 1

    return count


def approximate_pattern_match(Text, pattern, d):
    __d = []
    for i in range(0, len(Text) - len(pattern) + 1):
        if hammingdistance(Text[i:i + len(pattern)], pattern) <= d:
            __d.append(i)
    return __d


def frequent_words(Text, k):
    l = dict()
    for i in range(0, len(Text) - k + 1):
        kmer = Text[i:i + k]
        if kmer in l:
            l[kmer] += 1
        else:
            l[kmer] = 0
    return [kmer for kmer in l.keys() if l[kmer] == max(l.values())]


def frequent_words_with_mismatches_and_rcomplements(Text, k, d):
    __d = [0] * (4**k)
    for i in range(0, len(Text) - (k - 1)):
        pat = Text[i:i + k]
        close = neighbors(pat, d)
        for word in close:
            __d[patterntonumber(word)] += 1
        rpat = reversecomplement(pat)
        rclose = neighbors(rpat, d)
        for word in rclose:
            __d[patterntonumber(word)] += 1
    m = max(__d)
    values = [i for i, j in enumerate(__d) if j == m]
    return [numbertopattern(x, k) for x in values]


def motif_enumeration(Dna, k, d):
    k = int(k)
    d = int(d)
    sample = ""
    for word in Dna:
        sample += word
    close = [0] * 4**k
    patterns = []
    motifs = []
    for word in Dna:
        for i in range(0, (len(word) - k) + 1):
            neighborhood = neighbors(word[i:i + k], d)
            for i in range(0, len(neighborhood)):
                index = patterntonumber(neighborhood[i])
                close[index] = 1
        for i in range(0, 4**k):
            pat = numbertopattern(i, k)
            if close[i] == 1 and pat not in patterns:
                patterns.append(pat)
    for i in range(0, len(patterns)):
        count = 0
        for word in Dna:
            for j in range(0, (len(word) - k) + 1):
                if hammingdistance(word[j:j + k], patterns[i]) <= d:
                    count += 1
                    break
        if count == len(Dna):
            motifs.append(patterns[i])
    return motifs


def median_string(Dna, k):
    k = int(k)
    distance = float('inf')
    median = []
    for i in range(0, 4**k):
        d = sum([min([hammingdistance(numbertopattern(i, k), word[j:j + k])
                      for j in range(0, len(word) - k + 1)]) for word in Dna])
        if distance >= d:
            distance = d
            median.append([numbertopattern(i, k), distance])
    return [string for string in median
            if string[1] == min([m[1] for m in median])]


def most_probable_profile(text, profile, k):
    maxprob = -1
    for i in range(len(text) - k + 1):
        prob = 1
        for j, nucleotide in enumerate(text[i:i + k]):
            prob *= profile[nuc[nucleotide]][j]
        if prob > maxprob:
            maxprob = prob
            mostprob = text[i:i + k]
    return mostprob


def score(Motifs):
    prof = profile(Motifs)
    test = consensus(len(Motifs[0]), prof)

    return sum(
        [hammingdistance(test, Motifs[x]) for x in range(0, len(Motifs))])


def profile(Dna):
    t = len(Dna[0])
    __array = [[0.] * t, [0.] * t, [0.] * t, [0.] * t]
    for word in Dna:
        for index, key in enumerate(word):
            __array[nuc[key]][index] += 1
    dumby = []
    for i in range(0, 4):
        dumby.append([(x + 1) / (t + 1) for x in __array[i]])
    return dumby


def consensus(k, pr):
    motif = ""
    for i in range(0, k):
        column = [pr[base][i] for base, nucleotide in enumerate(nuc)]
        motif += sorted(list(nuc))[column.index(max(column))]
    return motif


def GreedyMotifSearch(Dna, k, t):
    k = int(k)
    bestscore = [t * k, None]
    for m in range(0, len(Dna[0]) - k + 1):
        motifs = [Dna[0][m:m + k]]
        profile_lst = profile(motifs)
        for n in range(1, t):
            motifs.append(most_probable_profile(Dna[n], profile, k))
            profile_lst = profile(motifs)
        currentscore = score(motifs)
        if currentscore < bestscore[0]:
            bestscore = [currentscore, motifs]

    return bestscore[1]


def RandomizedMotifSearch(Dna, k, t, i):
    from random import randint
    count = 0
    bestmotifs = [""] * t
    for j in range(0, t):
        r = randint(0, len(Dna[0]) - k)
        bestmotifs[j] = Dna[j][r:r + k]
    while count <= i:
        motifs = [""] * t
        for j in range(0, t):
            r = randint(0, len(Dna[0]) - k)
            motifs[j] = Dna[j][r:r + k]
        while True:
            profile_lst = profile(motifs)
            motifs = [most_probable_profile(Dna[x], profile, k)
                      for x in range(t)]
            if score(motifs) < score(bestmotifs):
                bestmotifs = motifs

            else:
                break

        count += 1
    return bestmotifs


def GibbsSampler(dna, k, t, N):
    from random import randint, randrange
    motifs = [dna[i][r:r + k]
              for i, r in enumerate([randint(0, len(dna[0]) - k)
                                     for a in range(t)])]
    bestmotifs = motifs
    for j in range(0, N):
        i = randrange(t)
        pr = profile(motifs[:i] + motifs[i + 1:])
        motifs[i] = most_probable_profile(dna[i], pr, k)
        if score(motifs) < score(bestmotifs):
            bestmotifs = motifs

    return bestmotifs


def Composition(Text, k):
    i = len(Text) - k + 1
    comp = [""] * i
    for n in range(0, i):
        comp[n] = Text[n:n + k]
    return sorted(comp)


def GenomePathString(genpath):
    k = len(genpath[0])
    string = ""
    string += genpath[0]
    for i in range(1, len(genpath)):
        j = 1
        while True:
            if genpath[i - 1][j:j + 2] == genpath[i][0:2]:
                break
            else:
                j += 1
        string += genpath[i][k - j:]
    return string


def OverlapGraph(kmers):
    k = len(kmers[0])
    n = k - 1
    adjlist = {}
    for kmer in kmers:
        if kmer[0:n] not in adjlist:
            adjlist[kmer[0:n]] = [kmer[1:k]]
        else:
            adjlist[kmer[0:n]].append(kmer[1:k])
    return adjlist


def deBruijnGraph(path, k):
    kmers = [path[i:i + k] for i in range(0, len(path) - k + 1)]
    return OverlapGraph(kmers)


def cycle(node, graph):
    from random import randint
    cycle = [node]
    while True:
        cycle.append(graph[node].pop(randint(0, len(graph[node]) - 1)))
        if len(graph[node]) == 0:
            del graph[node]
        if cycle[-1] in graph:
            node = cycle[-1]
        else:
            break
    return [cycle, graph]


def EulerianCycle(graph):
    from random import randint
    node = graph.keys()[randint(0, len(graph) - 1)]
    circuit, graph = cycle(node, graph)
    path = circuit
    while len(graph) > 0:
        for i in range(len(path)):
            if path[i] in graph:
                node = path[i]
                circuit, graph = cycle(node, graph)
                path = path[:i] + circuit + path[i + 1:]
    return path


def EulerianPath(graph):
    node = ""
    nodes = graph.keys()
    edges = [ll for l in graph.values() for ll in l]
    for n in nodes:
        if edges.count(n) < len(graph[n]):
            node = n
    if node == "":
        return EulerianCycle(graph)

    circuit, graph = cycle(node, graph)
    path = circuit
    while len(graph) > 0:
        for i in range(len(path)):
            if path[i] in graph:
                node = path[i]
                circuit, graph = cycle(node, graph)
                path = path[:i] + circuit + path[i + 1:]
    return path


def StringReconstruction(kmers):
    path = EulerianPath(OverlapGraph(kmers))
    if len(path[0]) == 1:
        return "".join(path)
    else:
        return path[0] + "".join([let[len(let) - 1] for let in path[1:]])


def StringSpelledByGappedPatterns(pairs, k, d):
    n = len(pairs)
    prefix, suffix = [StringReconstruction([pair[i] for pair in pairs])
                      for i in range(2)]
    checkpairs = [[prefix[i:i + k], suffix[i:i + k]]
                  for i in range(0, n + d - 1)]
    while prefix[k + d:] != suffix[:-(k + d)] or pairs != checkpairs:
        prefix, suffix = [StringReconstruction([pair[i] for pair in pairs])
                          for i in range(2)]
        checkpairs = [[prefix[i:i + k], suffix[i:i + k]]
                      for i in range(0, n + d - 1)]
    return prefix[:k + d] + suffix


def Contigs(kmers):
    from copy import deepcopy
    graph = OverlapGraph(kmers)
    nodes = [ll for l in [list(_[k]) for k in _ for _ in list(
        deepcopy(graph))] for ll in l]
    edges = [ll for l in graph.values() for ll in l]
    breaks = []
    for node in nodes:
        if len(graph[node]) != 1 or nodes.count(node) != edges.count(node):
            breaks.append(node)
    contigs = []
    for n in breaks:
        contig = [n]
        while True:
            n = graph[n].pop(0)
            contig.append(n)
            if n in breaks:
                break
            elif n not in nodes:
                break
        contigs.append(contig)
    return contigs


def kUniversalBinaryString(k):
    binary = int('0b' + '1' * k, 2)
    strings = []
    for i in range(0, binary + 1):
        num = str(bin(i))
        if len(num[2:]) < k:
            num = '0b' + '0' * (k - len(num[2:])) + num[2:]
            strings.append(num)
        else:
            strings.append(num)
    kmers = [num[2:] for num in strings]
    return EulerianPath(OverlapGraph(kmers))


def fragment(peptide, cyclic=False):
    fragments = []
    if cyclic == True:
        for string in [peptide[i:] + peptide[:i] for i in range(len(peptide))]:
            [fragments.append(string[0:j]) for j in range(1, len(peptide))]
        fragments.append(peptide)
    else:
        for i in range(len(peptide)):
            [fragments.append(peptide[i:i + j])
             for j in range(1, len(peptide) - i + 1)]
    return fragments


def theoreticalspectrum(peptide, cyclic=False):
    theory = [0]
    fragments = fragment(peptide, cyclic)
    [theory.append(sum([lookup.get_mass(ch) for ch in frag]))
     for frag in fragments]
    return sorted(theory)


def encoding_peptides(dna, peptide):
    rna = transcribe(dna)
    k = len(peptide) * 3
    rng = range(0, (len(dna) - k + (1 * len(peptide))) * 2)
    cond = lambda s: bool(translate(s) == peptide or translate(
        reversecomplementrna(s)) == peptide)
    return [reversetranscribe(rstrand)
            for rstrand in ["".join(rna)[i:i + k] for i in rng] if cond(rstrand)]


def peptide_sequencing(spectrum, cyclic=False):
    from itertools import permutations
    aminoacids = list(m for m in spectrum if m in lookup.acids)
    length = len(aminoacids)
    peptides = permutations(aminoacids)
    return [[lookup.get_mass(pep) for pep in peptide] for peptide in peptides if theoreticalspectrum(peptide, cyclic) == spectrum]


def spectrum_score(peptide, spectrum, cyclic=False):
    theoretical = theoreticalspectrum(peptide, cyclic)
    if spectrum[-1] < theoretical[-1]:
        return -1
    else:
        return sum([min(theoretical.count(protein), spectrum.count(protein))
                    for protein in set(spectrum)])


# need to translate this to just using masses
def leaderboard_sequencing(spectrum, N, cyclic=False):
    acids_masses = lookup.masses
    peptides = dict()
    best_seqs = list([0, set()])
    for a in acids_masses:
        score = spectrum_score(a, spectrum, cyclic)
        if score in peptides:
            peptides[score].add(a)
        else:
            peptides[score] = {a}
    while len(peptides) > 0:
        growths = dict()
        for vals in peptides.values():
            for seq in vals:
                for a in aa:
                    growth = seq + a
                    if mass(growth) in spectrum:
                        score = spectrum_score(growth, spectrum, cyclic)
                        if score != -1:
                            if score in growths:
                                growths[score].add(growth)
                            else:
                                growths[score] = {growth}

        if len(growths.keys()) > 0:
            peptides = growths
            max_score = max(peptides.keys())
        else:
            break

        if max_score >= best_seqs[0]:
            best_seqs[0] = max_score
            best_seqs[1] = peptides[max_score]
        total = sum([len(l) for l in peptides.values()])
        if sum([len(l) for l in peptides.values()]) > N:
            min_score = min(peptides.keys())
            while total > N:
                if len(peptides[min_score]) != 0:
                    peptides[min_score].pop()
                    total -= 1
                else:
                    del peptides[min_score]
                    min_score = min(peptides.keys())
    return '-'.join([str(acidsmasses[ch]) for ch in best_seqs[1].pop()])


def spectral_convolution_with_multiplicities(spectrum):
    convolution = dict()
    for i, peak in enumerate(spectrum):
        for mass in spectrum[:i]:
            diff = peak - mass
            if diff >= 57 and diff <= 200:
                if diff in convolution:
                    convolution[diff] += 1
                else:
                    convolution[diff] = 1
    return convolution


def spectral_convolution(spectrum):
    return list(
        filter(lambda x: x != None,
               [peak - mass if peak - mass > 0 else None
                for mass in spectrum[:i]
                for i, peak in enumerate(spectrum)])
    )


def convolution_cyclopeptide_sequencing(M, N, spectrum):
    convolution = spectral_convolution_with_multiplicities(spectrum)
    multiplicities = list()
    total = 0
    temp = sorted(convolution.values(), reversed=True)
    for i, multiplicity in enumerate(temp):
        total += multiplicity
        if total > M:
            multiplicities = set(temp[:i])
            break
    leaderboard_spectrum = leaderboard_sequencing(
        [k for k in convolution if convolution[k] in multiplicities], N, True)


if __name__ == "__main__":
    print(convolution_cyclopeptide_sequencing(20, 60, [int(
        x) for x in '57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493'.split(' ')]))
