def dp_change(m, d):
    mnc = [1 if i in d else 0 for i in range(0, m + 1)]
    for i in range(1, len(mnc)):
        if i not in d:
            mnc[i] = min(
                [mnc[i - c] + 1 if (i - c) > 0 else float('inf') for c in d])
    return mnc[m]


def _mtp_matrix(sink, down, right):  # manhattan tourist problem
    s = [[0 for i in range(0, sink[1] + 1)] for j in range(0, sink[0] + 1)]
    for i in range(1, len(s)):
        s[i][0] = s[i - 1][0] + down[i - 1][0]
    for j in range(1, len(s[0])):
        s[0][j] = s[0][j - 1] + right[0][j - 1]
    # to hell with being pythonic; listcomps are ALWAYS leagues faster and
    # still easy to read
    [s[i].__setitem__(j, max(
        [s[i - 1][j] + down[i - 1][j], s[i][j - 1] + right[i][j - 1]]))
     for j in range(1, sink[1] + 1) for i in range(1, sink[0] + 1)]
    return s


def mtp(sink, down, right):
    s = _mtp_matrix(sink, down, right)
    return s[sink[0]][sink[1]]


def lcs(v, w):
    s = [[0 for _ in range(len(w) + 1)] for __ in range(len(v) + 1)]
    for i in range(len(v)):
        for j in range(len(w)):
            if v[i] == w[j]:
                s[i + 1][j + 1] = s[i][j] + 1
            else:
                s[i + 1][j + 1] = max(s[i + 1][j], s[i][j + 1])
    lcs = ''
    i, j = len(v), len(w)
    while not any([i == 0, j == 0]):
        if s[i][j] == s[i - 1][j]:
            i -= 1
        elif s[i][j] == s[i][j - 1]:
            j -= 1
        else:
            lcs = v[i - 1] + lcs
            i -= 1
            j -= 1
    return lcs


def lcs_dag(source, sink, weighted_graph):
    seqs = []
    node = sink
    edges = [key for key in weighted_graph
             if node in list(weighted_graph[key].keys())]
    for edge in edges:
        score = weighted_graph[edge][node]
        seqs.append([score, edge, node])

    while not all([seq[1] == source for seq in seqs]):
        temp = []
        for seq in seqs:
            node = seq[1]
            edges = [key for key in weighted_graph
                     if node in list(weighted_graph[key].keys())]
            weights = [weighted_graph[edge][node] for edge in edges]
            if seq[1] == source:
                temp.append(seq)
                continue
            elif len(weights) == 0:
                continue
            else:
                for edge in edges:
                    score = weighted_graph[edge][node]
                    s = seq.copy()
                    s[0] += score
                    s.insert(1, edge)
                    temp.append(s)
        seqs = [item for item in temp]
    scores = [seq[0] for seq in seqs]
    return [seq for seq in seqs if seq[0] == max(scores)]
if __name__ == "__main__":
    print(lcs('GAGCAATT', 'ACTTAATT'))
