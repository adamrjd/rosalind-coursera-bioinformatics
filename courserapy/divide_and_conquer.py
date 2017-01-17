'''Written in Python36'''
from dataimport import Score


def show(matrix):
    for row in matrix:
        print('\t'.join([str(n) for n in row]))


def linear_scoring(v, w, sigma, score_dict):
    '''Linear scoring for first half'''
    # initialize and seed with penalties
    # this time we are actually splitting the vertical matrix
    # and adding one row at a time vertically
    s = [[-sigma * _] for _ in range(len(v) + 1)]
    backtrack = [None for _ in range(len(v) + 1)]

    # scoring
    for y in range(1, len(w) + 1):
        for x in range(len(v) + 1):
            if x == 0:
                s[x].append(-sigma * y)
            else:
                scores = (s[x - 1][1] - sigma,  # right
                          # compensate for len being 1 starting off
                          s[x][len(s[x]) - 1] - sigma,  # down
                          s[x - 1][0] + score_dict[v[x - 1]][w[y - 1]])  # diag

                score = max(scores)
                s[x].append(score)
                backtrack[x] = scores.index(score)

            if len(s[x]) > 2:
                s[x].pop(0)

    return s, backtrack


def middle_edge(v, w, sigma, score_dict):
    def get_edge(v, w, sigma, score_dict):
        '''
        The course actually sets it up so that the v string is on the vertical axis,
        and the w string is on the horizontal axis. I think this makes zero sense, and
        is much harder to reason about, so I flipped it.
        '''

        # will only ever catch len(w) == 2 b/c of outer return condition
        jcol = len(w) // 2 if len(w) > 2 else 1

        # scoring matrices for left and right
        rfoldr = lambda l: list(l)[::-1] if isinstance(l, list) else l
        s = linear_scoring(v, w[:jcol], sigma, score_dict)[0]
        r_s, backtrack = [
            [rfoldr(_) for _ in __] for __ in
            linear_scoring(v[::-1], w[jcol:][::-1], sigma, score_dict)]

        max_score = -float('inf')
        node = None
        for x in range(len(w) + 1):
            score = s[x][1] + r_s[x][0]
            if score > max_score:
                node, max_score = [(x, jcol), score]

        dirs = [tuple(map(sum, zip(node, _)))
                for _ in [(1, 0), (0, 1), (1, 1)]]

        return (node, dirs[backtrack[node[0]]])
    return (None, None) if len(w) <= 1 else get_edge(v, w, sigma, score_dict)


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
        def recurse(self, left, right, top, bottom):
            from dynamic_programming import global_alignment
            # right block empty, want to now align this segment
            if top == bottom:
                return [self.v[left:right], '-' * (right - left)]
            # left block empty, want to now align this segment
            elif left == right:
                return ['-' * (bottom - top), self.w[top:bottom]]

            # trivial cases exit recursion here
            elif bottom - top == 1 or right - left == 1:
                return global_alignment(self.sigma, self.score_dict,
                                        self.v[left:right], self.w[top:bottom])[1:]
            else:
                node, edge = middle_edge(
                    self.v[left:right], self.w[top:bottom], self.sigma, self.score_dict)

                # shift nodes from local source to global source
                node = tuple(map(sum, zip(node, [left, top])))
                edge = tuple(map(sum, zip(edge, [left, top])))

                # representation of middle node in v and w
                current = [
                    ['-', self.v[node[0]]][edge[1] - node[1]],
                    ['-', self.w[node[1]]][edge[0] - node[0]]
                ]

                # generate next subalignment from left and right subtrees
                A = recurse(self,
                            left, node[0], top, node[1])

                B = recurse(self,
                            edge[0], right,
                            edge[1], bottom)

                # zip left and right together recursively
                return [A[i] + current[i] + B[i] for i in range(2)]

        self.v, self.w = recurse(self, 0, len(self.v), 0, len(self.w))
        score = sum([-self.sigma if '-' in _ else self.score_dict[_[0]][_[1]]
                     for _ in zip(self.v, self.w)])

        return score, self.v, self.w

    def local_align(self):
        pass

if __name__ == "__main__":
    str1 = 'TWLNSACYGVNFRRLNPMNKTKWDCWTWVPMVMAAQYLCRIFIPVMDHWEFFGDWGLETWRLGIHDHVKIPNFRWSCELHIREHGHHFKTRFLKHNQFTQCYGLMPDPQFHRSYDVACQWEVTMSQGLMRFHRQNQIEKQRDRTSTYCMMTIGPGFTSNGYDPFVTITITPVQEPVENWFTPGGSMGFMIISRYMQMFFYLTRFSDMTYLVGVHCENYVCWNNVAKFLNGNLQGIFDQGERAYHQFVTWHSYSQYSRCSVGRYACEQAMSRVNSKMTWHWPIRDQGHEHFSEQYLSEKRNPPCNPRIGNAGQHFYEIHRIAHRVAMCNWAPQGQHPGGPTPHDVETCLWLWSLCLKGSDRGYVDRPWMFLADQLGEANLTLITMFHGCTRGCLMWFMDWEECVCSYSVVNPRCHGSEQWSVQNLGWRTCDTLISLWEPECDKHNTPPCLHWEFEDHPSQLRPVMMCDKYVQSIPTDAKWAWTYSKDFVISHWLIWTPIKLEECVFPQINRLWGTACNQGSQKIVIQNVWLRPSSFFQERSKCSDSSCILNVGGSNVNITGKETRTHVPILHMHEIDLISTASSGMRHNLILPHGMLMLHMNWHHSTRAMNPYSSLKLIPWTFQVCETDDRDQNVATHVADPCHKGEDQEIRCCKGGVDHQWKGDRMWMMCMPDMNYVKQDQAPSGTCEGACENYPADKDKCYMIFTIVFDYRRCTKKVCIWISGFPVDAFNLISIANAGFFCCWLEPTELKWRRTFYLGKGTQGWMCTFPHRNIIPVIICAGFGRWVQGEVPFRPVAQISAHSSDRRQGHHPPGTNMCHDYGDQYPIKRVGMQVEEDDGASYCDCAADWKLADMYEADHLSIGVIDFTDWIYPKNGGIWSEIIKSHFHWYHWETPQNTVGAFNTIVGINGSDMCIYHGNTQWEFGWCWKWLNHGHMRNQGPCHLGILEGRISKFAQVTSWWWQTKHDKDWSIEPYGRHWGEAGRPYTYNYCWMRWAIVYNHGNVISVELVPFMDEYPGKCNKEDVQFELFSPMQA'
    str2 = 'LWFKFLQCIFQYFKDQQETNCIWTFSPFSEHICQRVCQVYWNWNTPSSRTSDPRELFANSTIHNNRCGEWRYMFYHTRTLVQTAPLMKETLHSDGKHSMYCEQRHFFRSSYLIKVNYDVSHYLELYTFSEIPWKLTTHGWDGFSWFLLVNSCCTFDIDGKCGILSQCGMSRAFRTRQEDAYHFQTSLMHLHLHLHVQEGKHEKADLFAQFYNMLPMHGGTCGRNTEPSDLFDSATMNKYMAEHPASCKACPNVSKECFVYWWSHDFTKKHKLIEFSCGRDTGQTTQRTWNVDENEGGKWIWRFHYFMRAKALQIDPKFKPYWNEPRAIMRPGHVTAAPCICAQHSQNETAVCNRDQMHIHAIEFQQYHSRAFGEVQTWCDIGKENENDFIYEQHWWLVGGTEGMAGVIWKFVCARCRTQDCDFWKTCLTYSAQPMMKVYDTIFYVNSINPWEFEDHPSQCDKCVQSIPTDAKYAICGKFVISHWLYWTPQKFEECVHNNVRCAPMGNRLWGTACMVIQNVWLRPSMGSHFSCILNVGGSNINIQGKETWTHVPILHMHEIDLISTASSGMETCKPCFLSGPTIHMGFSYEIRAQPYSRDYFCMDWMQEADEVDHNRCETVQPTLPLLQQFEWKTSCMGQRWITIFCDHCQIVCFSTFFCVMPTFLPNTSILDKFYCIYLSISWTHYCNVHALGFIMRLHYSYMGWKEHKRMHAWDIGLDELWAQEGIQRAQLWCGDEFEVAKYPEWITEARTAIATRPWFHNCYIKPWWIREKHLWFGKESKLDHGHRGAMFTPVANDNTEWMHHWYMFCWAGSKNRLKRQIKEKLIFIIKFMITEFGLFLMIDYTQCYIAWMWAYTGIACYIDWEKCLKHDLTTTDLGCCVYRLFKWYEVRHRAPPQVNTRLPWSQIPMVAIQCNIVDECKEQWHFSYKASFVVEYLCPGCCTNGNRWQWYQVKETPFMYAFAASIFGFHHENLVVFITGSVTIPNGLFGCIAWTSPKPVQKTPASANTIIAYDKCILMG'
    penalty = 5
    score_matrix = Score().BLOSUM62
    print(middle_edge(str1, str2, penalty, score_matrix))
    '''answer = LinearSpaceAlignment(
        str1, str2, penalty, score_matrix).global_align()
    print(answer)'''
