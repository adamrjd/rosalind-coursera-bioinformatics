'''really need to refactor this code; could do way better'''

from scripts import *
nuc = {nucleotide: index for index, nucleotide in enumerate('ACGT')}


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


def consensus(k, profile):
    motif = ""
    for i in range(0, k):
        column = [profile[base][i] for base, nucleotide in enumerate(nuc)]
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
    count = 0
    bestmotifs = [""] * t
    for j in range(0, t):
        r = random.randint(0, len(Dna[0]) - k)
        bestmotifs[j] = Dna[j][r:r + k]
    while count <= i:
        motifs = [""] * t
        for j in range(0, t):
            r = random.randint(0, len(Dna[0]) - k)
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
    motifs = [dna[i][r:r + k]
              for i, r in enumerate([random.randint(0, len(dna[0]) - k)
                                     for a in range(t)])]
    bestmotifs = motifs
    for j in range(0, N):
        i = random.randrange(t)
        profile = profile(motifs[:i] + motifs[i + 1:])
        motifs[i] = most_probable_profile(dna[i], profile, k)
        if Score(motifs) < Score(bestmotifs):
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
    cycle = [node]
    while True:
        cycle.append(graph[node].pop(random.randint(0, len(graph[node]) - 1)))
        if len(graph[node]) == 0:
            del graph[node]
        if cycle[-1] in graph:
            node = cycle[-1]
        else:
            break
    return [cycle, graph]


def EulerianCycle(graph):
    node = graph.keys()[random.randint(0, len(graph) - 1)]
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
    edges = flatten(graph.values())
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
    graph = OverlapGraph(kmers)
    dumby = OverlapGraph(kmers)
    nodes = flatten([[key for dumby[key] in dumby[key]] for key in dumby])
    edges = flatten(graph.values())
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

if __name__ == "__main__":
    with open('text.txt', 'r') as f:
        read_data = [l.strip('\n()') for l in f.readlines()]
    pairs = [_.split('|') for _ in read_data]
    print(StringSpelledByGappedPatterns(pairs, 3, 1))
