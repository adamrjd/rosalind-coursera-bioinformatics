import os


class BLOSUM62(object):
    """The BLOSUM62 scoring matrix class."""

    def __init__(self):
        """Initialize the scoring matrix."""
        with open(os.path.join(os.path.dirname(__file__), 'datfiles/BLOSUM62_dumby.txt')) as input_data:
            items = [line.strip().split() for line in input_data.readlines()]
            self.scoring_matrix = {
                (item[0], item[1]): int(item[2]) for item in items}

    def __getitem__(self, pair):
        """Returns the score of the given pair of protein."""
        return self.scoring_matrix[pair[0], pair[1]]


def middle_column_score(v, w, scoring_matrix, sigma):
    '''Returns the score of the middle column for the alignment of v and w.'''
    # Initialize the score columns.
    S = [[i*j*sigma for j in xrange(-1, 1)] for i in xrange(len(v)+1)]
    S[0][1] = -sigma
    backtrack = [0]*(len(v)+1)

    # Fill in the Score and Backtrack matrices.
    for j in xrange(1, len(w)/2+1):
        for i in xrange(0, len(v)+1):
            if i == 0:
                S[i][1] = -j*sigma
            else:
                scores = [S[i-1][0] + scoring_matrix[v[i-1], w[j-1]], S[i][0] - sigma, S[i-1][1] - sigma]
                S[i][1] = max(scores)
                backtrack[i] = scores.index(S[i][1])

        if j != len(w)/2:
            S = [[row[1]]*2 for row in S]

    return [row[1] for row in S], backtrack


def middle_edge(v, w, scoring_matrix, sigma):
    '''Returns the middle edge in the alignment graph of v and w.'''

    # Get the score of the middle column from the source to the middle.  The backtrack matrix is unnecessary here.
    source_to_middle = middle_column_score(v, w, scoring_matrix, sigma)[0]

    # Get the score of the middle column from the middle to sink.  Reverse the order as the computations are done in the opposite orientation.
    middle_to_sink, backtrack = map(lambda l: l[::-1], middle_column_score(v[::-1], w[::-1]+['', '$'][len(w) % 2 == 1 and len(w) > 1], scoring_matrix, sigma))

    # Get the componentwise sum of the middle column scores.
    scores = map(sum, zip(source_to_middle, middle_to_sink))
    with open('output.txt','w') as f:
        for score in scores:
            f.write(str(score)+'\n')

    # Get the position of the maximum score and the next node.
    max_middle = max(xrange(len(scores)), key=lambda i: scores[i])

    if max_middle == len(scores) - 1:
        next_node = (max_middle, len(w)/2 + 1)
    else:
        next_node = [(max_middle + 1, len(w)/2 + 1), (max_middle, len(w)/2 + 1), (max_middle + 1, len(w)/2),][backtrack[max_middle]]

    return (max_middle, len(w)/2), next_node


if __name__ == '__main__':
    v,w='PLEASNTLY','MEASNLY'
    middle = middle_edge(v, w, BLOSUM62(), 5)
    print(middle)