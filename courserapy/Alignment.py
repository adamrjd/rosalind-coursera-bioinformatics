'''
Module written for Coursera Bioinformatics III
Written by: Adam Dornford
Written in: Python35
'''

from Parsing import Score


def show(matrix):
    for row in matrix:
        print('\t'.join([str(n) for n in row]))


class DynamicAlignment(object):

    def __init__(self):
        pass

    @classmethod
    def scoring(cls, v, w, sigma, score_dict, local=False):
        i, j = [len(_) for _ in [v, w]]
        if local:
            max_pos = (0, [0, 0])
        score = lambda x, y: score_dict[v[x - 1]][w[y - 1]]
        # assemble matrices and account for indels
        s = [[0 for _ in range(j + 1)] for __ in range(i + 1)]
        backtrack = [[0 for _ in range(j + 1)] for __ in range(i + 1)]
        # fill in penalties
        [s[x].__setitem__(0, -sigma * x) for x in range(1, i + 1)]
        [s[0].__setitem__(y, -sigma * y) for y in range(1, j + 1)]
        # create alignment scoring matrix
        for x in range(1, i + 1):
            for y in range(1, j + 1):
                l = [
                    s[x - 1][y] - sigma,  # right
                    s[x][y - 1] - sigma,  # down
                    s[x - 1][y - 1] + score(x, y)]  # diag
                s[x][y], backtrack[x][y] = max(l), l.index(max(l))
                if local:
                    if s[x][y] > max_pos[0]:
                        max_pos[0] = s[x][y]
                        max_pos[1] = [x, y]

        if local:
            return s, backtrack, max_pos
        return s, backtrack

    @classmethod
    def global_alignment(cls, sigma, score_dict, v, w):
        i, j = [len(v), len(w)]
        score = 0
        s, backtrack = DynamicAlignment.scoring(sigma, score_dict, v, w)
        indel = lambda word, i: word[:i] + '-' + word[i:]
        score = s[i][j]
        while i * j != 0:
            if backtrack[i][j] == 0:
                i -= 1
                w = indel(w, j)
            elif backtrack[i][j] == 1:
                j -= 1
                v = indel(v, i)
            else:
                i -= 1
                j -= 1
        for _ in range(i):
            w = indel(w, 0)
        for _ in range(j):
            v = indel(v, 0)
        return score, v, w

    '''def local_alignment(sigma, score_dict, v, w):
        # this is broken
        i, j = [len(v), len(w)]
        a = [0, v, w]
        s, backtrack, max_pos = scoring(sigma, score_dict, v, w, True)
        a[0] = max_pos[0]
        i, j, a[1], a[2] = align(backtrack, v, w, max_pos[1])
        return a'''

    @classmethod
    def edit_distance(cls, v, w, sigma, score_dict):
        i, j = [len(v), len(w)]

        # assemble matrices and account for indels
        s = [[0 for _ in range(j + 1)] for __ in range(i + 1)]
        for x in range(1, i + 1):
            s[x][0] = s[x - 1][0] + sigma
        for y in range(1, j + 1):
            s[0][y] = s[0][y - 1] + sigma

        for x in range(1, i + 1):
            for y in range(1, j + 1):
                l = [s[x - 1][y] + sigma, s[x][y - 1] + sigma,
                     s[x - 1][y - 1] + score_dict[v[x - 1]][w[y - 1]]]
                s[x][y] = min(l)
        return s[i][j]

    @classmethod
    def affine_gap_alignment(cls, v, w, sigma, epsilon, score_dict):
        i, j = [len(v), len(w)]
        score = 0
        # lambdas
        scoreij = lambda x, y: score_dict[v[x - 1]][w[y - 1]]
        indel = lambda seq, i: seq[:i] + '-' + seq[i:]
        lowerij = lambda x, y: [s[0][y - 1][x] - epsilon,
                                s[1][y - 1][x] - sigma, -float('inf')]
        middleij = lambda x, y: [s[0][y][x], s[1][
            y - 1][x - 1] + scoreij(x, y), s[2][y][x]]
        upperij = lambda x, y: [-float('inf'), s[1][y]
                                [x - 1] - sigma, s[2][y][x - 1] - epsilon]

        # initialize upper, middle, and lower backtrack matrices
        b = [[[0 for _ in range(i + 1)] for __ in range(j + 1)]
             for ___ in range(3)]
        s = [[[0 for _ in range(i + 1)] for __ in range(j + 1)]
             for ___ in range(3)]

        # seed with penalties
        [[s[0][0].__setitem__(x, -10 * sigma), s[1][0].__setitem__(x, -sigma - (x - 1) * epsilon),
          s[2][0].__setitem__(x, -sigma - (x - 1) * epsilon)] for x in range(1, i + 1)]
        [[s[0][y].__setitem__(0, -sigma - (y - 1) * epsilon), s[1][y].__setitem__(0, -sigma - (
            y - 1) * epsilon), s[2][y].__setitem__(0, -10 * sigma)] for y in range(1, j + 1)]

        # execute scoring and backtracking
        for x in range(1, i + 1):
            for y in range(1, j + 1):
                vals = lowerij(x, y)
                s[0][y][x], b[0][y][x] = [max(vals), vals.index(max(vals))]
                vals = upperij(x, y)
                s[2][y][x], b[2][y][x] = [max(vals), vals.index(max(vals))]
                # upper and lower before middle
                vals = middleij(x, y)
                s[1][y][x], b[1][y][x] = [max(vals), vals.index(max(vals))]

        point = 1
        score = s[1][j][i]
        # dynamically trace back
        while i * j != 0:
            if point == 0:  # lower traceback
                point = b[0][j][i]
                v = indel(v, i)
                j -= 1
            elif point == 1:  # middle traceback
                if b[1][j][i] == 0:
                    point = b[1][j][i]
                    j -= 1
                    v = indel(v, i)
                if b[1][j][i] == 2:
                    point = b[1][j][i]
                    i -= 1
                    w = indel(w, j)
                else:
                    i -= 1
                    j -= 1
                    point = b[1][j][i]
            elif point == 2:  # upper traceback
                point = b[2][j][i]
                w = indel(w, j)
                i -= 1
        # trace back to origin
        for _ in range(i):
            w = indel(w, 0)
        for _ in range(j):
            v = indel(v, 0)

        return score, v, w

    @classmethod
    def multiple_sequence_alignment(cls, v, w, u):
        from math import ceil
        l, m, n = [len(_) for _ in (v, w, u)]
        align_score = 0
        # lambdas
        indel = lambda seq, p: seq[:p] + '-' + seq[p:]

        # initialize matrices
        s = [[[0 for _ in range(n + 1)] for __ in range(m + 1)]
             for ___ in range(l + 1)]
        backtrack = [[[0 for _ in range(n + 1)]
                      for __ in range(m + 1)] for ___ in range(l + 1)]
        for x in range(1, l + 1):
            for y in range(1, m + 1):
                for z in range(1, n + 1):
                    scores = [
                        s[x - 1][y - 1][z - 1] +
                        int(v[x - 1] == w[y - 1] == u[z - 1]),
                        s[x][y - 1][z - 1],
                        s[x - 1][y][z - 1],
                        s[x - 1][y - 1][z],
                        s[x][y][z - 1],
                        s[x][y - 1][z],
                        s[x - 1][y][z],
                    ]
                    s[x][y][z] = max(scores)
                    backtrack[x][y][z] = scores.index(s[x][y][z])
        align_score = ceil(s[l][m][n])

        # backtracking
        while l * m * n != 0:
            trace = backtrack[l][m][n]
            if trace == 6:
                l -= 1
                w = indel(w, m)
                u = indel(u, n)
            elif trace == 5:
                v = indel(v, l)
                m -= 1
                u = indel(u, n)
            elif trace == 4:
                v = indel(v, l)
                w = indel(w, m)
                n -= 1
            elif trace == 3:
                l -= 1
                m -= 1
                u = indel(u, n)
            elif trace == 2:
                l -= 1
                w = indel(w, m)
                n -= 1
            elif trace == 1:
                v = indel(v, l)
                m -= 1
                n -= 1
            elif trace == 0:
                l -= 1
                m -= 1
                n -= 1

        # finish insertions for all 3 strings if necessary
        for _ in range(l):
            v = indel(v, 0)
        for _ in range(m):
            w = indel(w, 0)
        for _ in range(n):
            u = indel

        return align_score, v, w, u


class LinearAlignment(object):
    '''
    Because of the way this problem is formulated, linear space alignment imo is best handled as a
    class...that way, there are no side effects from mutating v and w during recursion (at least
    in py3; in py2 you could have used the leaky 'lexical' scope bindings to simulate this same
    concept)
    '''

    def __init__(self, v, w, sigma, score_dict):
        self.v = v
        self.w = w
        self.sigma = sigma
        self.score_dict = score_dict

    @classmethod
    def linear_scoring(cls, v, w, sigma, score_dict):
        '''Creates '''
        # initialize and seed with penalties
        # dictionaries are much more efficient here since we
        # will still end up with very large data structures
        s = {(_, 0): -sigma * _ for _ in range(len(v) + 1)}
        backtrack = {_: 1 for _ in range(len(v) + 1)}

        # scoring
        for y in range(1, len(w) + 1):
            for x in range(len(v) + 1):
                if x == 0:
                    s[(x, y)] = -sigma * y
                else:
                    scores = (
                        s[(x - 1, y - 1)] +
                        score_dict[v[x - 1]][w[y - 1]],  # diag
                        s[(x, y - 1)] - sigma,  # right
                        s[(x - 1, y)] - sigma)  # down

                    score = max(scores)
                    s[(x, y)] = score
                    backtrack[x] = scores.index(score)
                if y - 2 >= 0 and len(v) > 2:
                    s.pop((x, y - 2))

        return s, backtrack

    @classmethod
    def middle_edge(cls, v, w, sigma, score_dict):
        def get_edge(v, w, sigma, score_dict):
            jcol = len(w) // 2

            # scoring matrices for left and right
            s = LinearAlignment.linear_scoring(
                v, w[:jcol], sigma, score_dict)[0]
            rev_s, rev_backtrack = LinearAlignment.linear_scoring(
                v[::-1], w[-(len(w) - jcol):][::-1], sigma, score_dict)

            max_score = -float('inf')
            node = None
            # len(v) // 2 + 1 approximate bounds for 'longest path' since sometimes
            # the max score occurs after the node with the longest path
            # also, need to account for special case where v very short and w
            # long

            for x in range(len(v) // 2 + 1):
                score = s[(x, jcol)] + rev_s[(len(v) - x, jcol)]
                if score > max_score:
                    node, max_score = [(x, jcol), score]

            # diagonal, right, down directions for edge from node
            dirs = [tuple(map(sum, zip(node, _)))
                    for _ in [(1, 1), (1, 0), (0, 1)]]

            # get edge from backtrack matrix
            edge = dirs[rev_backtrack[len(v) - node[0]]]
            return (node, edge)
        return (None, None) if len(v) == len(w) == 1 else tuple(
            tuple(reversed(_))  # reverse v and w if v too small and w large
            for _ in get_edge(w, v, sigma, score_dict)
        ) if len(v) <= 2 and len(w) > 1 else get_edge(v, w, sigma, score_dict)

    def global_align(self):
        def recurse(self, positions, history=None):
            from Utils import Tree
            left, right, top, bottom = positions
            if history is None:
                history = Tree()

            # base cases
            if top == bottom:
                return
            if left == right:
                return

            # get the current middle edge and shift with currnet source
            node, edge = [
                tuple(map(sum, zip(_, [left, top])))
                if _ is not None else None
                for _ in LinearAlignment.middle_edge(
                    self.v[left:right],
                    self.w[top:bottom],
                    self.sigma,
                    self.score_dict)]

            # if there is no middle edge, find last edge and shift all following
            # edges if there was indel
            if (node, edge) == (None, None):
                node, edge = (left, top), (right, bottom)
            else:
                history.extend([node, edge])

            # check for indel; shift edges if indel
            for _, __ in enumerate([edge[_] - node[_] for _ in range(1, -1, -1)]):
                if __ == 0:
                    history.pop(edge).extend(item[_] - 1 for item in history)

            recurse(self, [left, node[0], top, node[1]], history, node)

            recurse(self, [edge[0], right, edge[1], bottom], history, edge)

            return history

        # entry point for recursion
        lst = list(recurse(self, [0, len(self.v), 0, len(self.w)]))

        # initialize alignment
        A = ''
        B = ''
        align = lambda t1, t2: [
            ['-', self.v[t1[0]]][t2[1] - t1[0]], ['-', self.w[t1[1]]][t2[0] - t1[0]]]

        # from list of edges, get alignment
        while len(lst) > 0:
            a, b, *lst = lst
            temp = align(a, b)
            A += temp[0]
            B += temp[1]

        # set alignment
        self.v, self.w = A, B

        # score alignment
        score = sum([-self.sigma if '-' in _ else self.score_dict[_[0]][_[1]]
                     for _ in zip(self.v, self.w)])

        return score, self.v, self.w

    def local_align(self):
        '''to do'''
        pass

if __name__ == "__main__":
    with open('input.txt', 'r') as f:
        str1, str2 = [line.strip() for line in f.readlines()]
    penalty = 5
    score_matrix = Score().BLOSUM62
    with open('output.txt', 'w') as f:
        for _ in LinearAlignment(
                str1, str2, penalty, score_matrix).global_align():
            f.write(str(_) + '\n')
