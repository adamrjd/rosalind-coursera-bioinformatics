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
        r_s, r_backtrack = linear_scoring(
            v[::-1], w[-(len(w) - jcol):][::-1], sigma, score_dict)

        max_score = -float('inf')
        node = None
        # len(v) // 2 + 1 approximate bounds for 'longest path' since sometimes
        # the max score occurs after the node with the longest path
        # also, need to account for special case where v very short and w long

        for x in range(len(v) // 2 + 1):
            score = s[(x, jcol)] + r_s[(len(v) - x, jcol)]
            if score > max_score:
                node, max_score = [(x, jcol), score]

        # diagonal, right, down directions for edge from node
        dirs = [tuple(map(sum, zip(node, _)))
                for _ in [(1, 1), (1, 0), (0, 1)]]

        # get edge from backtrack matrix
        edge = dirs[r_backtrack[len(v) - node[0]]]
        return (node, edge)
    return ((0, 0), (1, 1)) if len(v) == len(w) == 1 else tuple(
        tuple(reversed(_))
        for _ in get_edge(w, v, sigma, score_dict)
    ) if len(v) <= 2 and len(w) > 1 else get_edge(v, w, sigma, score_dict)


class LinearSpaceAlignment(object):
    '''
    Because of the way this problem is formulated, linear space alignment imo is best handled as a
    class...that way, there are no side effects from mutating v and w during recursion (at least
    in py3; in py2 you could have used the leaky 'lexical' scope bindings to simulate this same
    concept)

    This works perfectly, but because the middle edge problem is formulated poorly middle_edge
    breaks for edge cases, e.g. 'PL','M' aligns to 'PL','-M' but there is no "middle edge";
    'L', 'M' should be 'L', '-M' for a global alignment of a particular string but the subalignment
    is of course 'L','M'; 'LY','Y' aligns to 'LY','-Y' but again only two edges not three
    '''

    def __init__(self, v, w, sigma, score_dict):
        self.v = v
        self.w = w
        self.sigma = sigma
        self.score_dict = score_dict

    def global_align(self):
        def recurse(self, left, right, top, bottom):
            # right block empty, want to now align this segment
            if top == bottom:
                subalignment = [self.v[left:right], '-' * (right - left)]
                return subalignment
            # left block empty, want to now align this segment
            elif left == right:
                subalignment = ['-' * (bottom - top), self.w[top:bottom]]
                return subalignment

            else:
                # get middle node and edge and shift with current source
                node, edge = [
                    tuple(map(sum, zip(_, [left, top])))
                    if _ is not None else None
                    for _ in middle_edge(
                        self.v[left:right],
                        self.w[top:bottom],
                        self.sigma,
                        self.score_dict)]

                # representation of middle node in v and w
                subalignment = [
                    ['-', self.v[node[0]]][edge[1] - node[1]],
                    ['-', self.w[node[1]]][edge[0] - node[0]]
                ]

                # generate next subalignment from left and right subtrees
                A = recurse(self,
                            left, node[0], top, node[1])

                # shift edge with current alignment if there is an indel
                if '-' in subalignment:
                    edge = tuple(edge[_] - 1 if subalignment[_] == '-' else edge[_] + 1
                                 for _ in range(2))

                B = recurse(self,
                            edge[0], right,
                            edge[1], bottom)

                # create alignment from subalignments
                return [A[i] + subalignment[i] + B[i] for i in range(2)]

        self.v, self.w = recurse(self, 0, len(self.v), 0, len(self.w))
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
