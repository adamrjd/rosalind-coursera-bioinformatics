'''Written in Py35'''
from dataimport import Score


def show(matrix):
    for row in matrix:
        print('\t'.join([str(n) for n in row]))


def linear_scoring(v, w, sigma, score_dict):
    '''Linear scoring for first half'''
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
                    s[(x - 1, y - 1)] + score_dict[v[x - 1]][w[y - 1]],  # diag
                    s[(x, y - 1)] - sigma,  # right
                    s[(x - 1, y)] - sigma)  # down

                score = max(scores)
                s[(x, y)] = score
                backtrack[x] = scores.index(score)
            if y - 2 >= 0 and len(v) > 2:
                s.pop((x, y - 2))

    return s, backtrack


def middle_edge(v, w, sigma, score_dict):
    def get_edge(v, w, sigma, score_dict):
        jcol = len(w) // 2

        # scoring matrices for left and right
        s = linear_scoring(v, w[:jcol], sigma, score_dict)[0]
        rev_s, rev_backtrack = linear_scoring(
            v[::-1], w[-(len(w) - jcol):][::-1], sigma, score_dict)

        max_score = -float('inf')
        node = None
        # len(v) // 2 + 1 approximate bounds for 'longest path' since sometimes
        # the max score occurs after the node with the longest path
        # also, need to account for special case where v very short and w long

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


class LinearSpaceAlignment(object):
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

    def global_align(self):
        def recurse(self, positions, history=None, parent_node=None):
            from utils import Tree
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
                for _ in middle_edge(
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
                    history.pop([edge]).extend(item[i] - 1 for item in history)

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
        for _ in LinearSpaceAlignment(
                str1, str2, penalty, score_matrix).global_align():
            f.write(str(_) + '\n')
