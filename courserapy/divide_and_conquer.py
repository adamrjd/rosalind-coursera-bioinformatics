'''Written in Python36'''
from dataimport import Score


def show(matrix):
    for row in matrix:
        print('\t'.join([str(n) for n in row]))


def linear_scoring(v, w, sigma, score_dict):
    '''Linear scoring for first half'''
    i, j = len(v), len(w)
    # initialize and seed with penalties
    s = [[] for _ in range(j + 1)]
    # scoring
    for x in range(i + 1):
        for y in range(j + 1):
            if y == 0 or x == 0:
                s[y] = [-sigma * y, - sigma * x]
            else:
                if v[x - 1] == w[y - 1]:
                    score = s[y - 1][0] + score_dict[v[x - 1]][w[y - 1]]
                else:
                    score = max(s[y][0], s[y - 1][1])
                s[y].append(score)
                s[y].pop(0)
    return s


def middle_edge(v, w, sigma, score_dict):
    def get_edge(v, w, sigma, score_dict):
        v, w = w, v
        i, j = len(v) + 1, len(w) + 1
        jcol = int((i - 1) / 2)

        # scoring matrices for left and right
        s = linear_scoring(v[:jcol], w, sigma, score_dict)
        r_s = [_ for _ in reversed(linear_scoring(
            v[jcol:][::-1], w[::-1], sigma, score_dict))]
        max_score = -float('inf')
        # right, diag, down
        right_nodes = [float('inf'), float('inf'), float('inf')]
        for y in range(jcol + 1):
            node = s[y][1] + r_s[y][0]
            if node > max_score:
                l_node = (y, jcol)  # since we reversed v and w
                max_score = node
                right_nodes[0] = r_s[y][1]
                if j > y + 1:
                    right_nodes[1] = r_s[y + 1][1]
                    right_nodes[2] = s[y + 1][1] + r_s[y + 1][0]
                else:
                    right_nodes[1], right_nodes[2] = [float('inf')] * 2
        if right_nodes[1] < right_nodes[2]:  # diag
            r_node = (l_node[0] + 1, l_node[1] + 1)
        else:  # right
            r_node = (l_node[0] + 1, l_node[1])
        return [l_node, r_node]
    return None if len(w) <= 1 else get_edge(v, w, sigma, score_dict)


class LinearSpaceAlignment(object):
    '''
    Because of the way this problem is formulated, linear space alignment imo is best handled as a class...
    that way, there are no side effects from mutating v and w during recursion
    '''

    def __init__(self, v, w, sigma, score_dict):
        self.v = v
        self.w = w
        self.sigma = sigma
        self.score_dict = score_dict

    def global_align(self):
        def recurse(self, top, bottom, left, right):
            from dynamic_programming import global_alignment
            # right block empty, want to now align this segment
            if left == right:
                return [self.v[top:bottom], '-' * (bottom - top)]
            # left block empty, want to now align this segment
            elif top == bottom:
                return ['-' * (right - left), self.w[left:right]]

            # trivial cases
            elif bottom - top == 1 or right - left == 1:
                return global_alignment(self.sigma, self.score_dict,
                                        self.v[top:bottom], self.w[left:right])[1:]
            else:
                mid_node, next_node = middle_edge(
                    self.v[top:bottom], self.w[left:right], self.sigma, self.score_dict)

                # shift nodes
                mid_node = tuple(map(sum, zip(mid_node, [top, left])))
                next_node = tuple(map(sum, zip(next_node, [top, left])))

                # representation of middle node in v and w
                current = [['-', self.v[mid_node[0] % len(self.v)]][next_node[0] - mid_node[0]], [
                    '-', self.w[mid_node[1] % len(self.w)]][next_node[1] - mid_node[1]]]

                # generate alignment from left and right subtrees
                A = recurse(self,
                            top, mid_node[0], left, mid_node[1])
                B = recurse(self,
                            next_node[0], bottom, next_node[1], right)
                return [A[i] + current[i] + B[i] for i in range(2)]

        self.v, self.w = recurse(self, 0, len(self.v), 0, len(self.w))
        score = sum([-self.sigma if '-' in _ else self.score_dict[_[0]][_[1]]
                     for _ in zip(self.v, self.w)])

        return score, self.v, self.w

    def local_align(self):
        pass

if __name__ == "__main__":
    str2 = 'PLEASANTLY'
    str1 = 'MEANLY'
    penalty = 5
    score_matrix = Score().BLOSUM62
    answer = LinearSpaceAlignment(
        str1, str2, penalty, score_matrix).global_align()

    print(answer[1:][::-1], answer[0])
