'''
Module written for Coursera Bioinformatics III
Written by: Adam Dornford
Written in: Python35
'''


def scoring(sigma, score_dict, v, w, local=False):
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
            l = [s[x - 1][y] - sigma, s[x][y - 1] -
                 sigma, s[x - 1][y - 1] + score(x, y)]
            s[x][y], backtrack[x][y] = max(l), l.index(max(l))
            if local:
                if s[x][y] > max_pos[0]:
                    max_pos[0] = s[x][y]
                    max_pos[1] = [x, y]

    if local:
        return s, backtrack, max_pos
    return s, backtrack


def global_alignment(sigma, score_dict, v, w):
    i, j = [len(v), len(w)]
    score = 0
    s, backtrack = scoring(sigma, score_dict, v, w)
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
    #this is broken
    i, j = [len(v), len(w)]
    a = [0, v, w]
    s, backtrack, max_pos = scoring(sigma, score_dict, v, w, True)
    a[0] = max_pos[0]
    i, j, a[1], a[2] = align(backtrack, v, w, max_pos[1])
    return a'''


def edit_distance(v, w, sigma, score_dict):
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


def affine_gap_alignment(sigma, epsilon, score_dict, v, w):
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


def multiple_sequence_alignment(v, w, u):
    l, m, n = [len(_) for _ in (v, w, u)]
    align_score = 0
    # lambdas
    scores_ijk = lambda i, j, k: [[
        s[i][j][k - 1],
        s[i][j - 1][k],
        s[i - 1][j][k],
        s[i][j - 1][k - 1],
        s[i - 1][j][k - 1],
        s[i - 1][j - 1][k],
        s[i - 1][j - 1][k - 1] + 1 if v[k - 1] == w[j - 1] == u[i - 1] else 0]]
    indel = lambda seq, p: seq[:p] + '-' + seq[p:]

    # initialize scoring matrix
    s = [[[0 for _ in range(l + 1)] for __ in range(m + 1)]
         for ___ in range(n + 1)]
    backtrack = [[[0 for _ in range(l + 1)]
                  for __ in range(m + 1)] for ___ in range(n + 1)]
    for x in range(1, l + 1):
        for y in range(1, m + 1):
            for z in range(1, n + 1):
                scores = scores_ijk(z, y, x)
                s[z][y][x] = max(scores)
                backtrack[z][y][x] = scores.index(s[z][y][x])
    align_score = s[n][m][l]

    # backtracking
    while l * m * n != 0:
        dex = backtrack[n][m][l]
        if dex == 2:
            w = indel(w, m)
            u = indel(u, n)
            l -= 1
        elif dex == 1:
            v = indel(v, l)
            u = indel(u, n)
            m -= 1
        elif dex == 0:
            v = indel(v, l)
            w = indel(w, m)
            n -= 1
        elif dex == 5:
            u = indel(u, n)
            l -= 1
            m -= 1
        elif dex == 4:
            w = indel(w, m)
            l -= 1
            n -= 1
        elif dex == 3:
            v = indel(v, l)
            m -= 1
            n -= 1
        elif dex == 6:
            l -= 1
            m -= 1
            n -= 1
    # backtrack to 0
    for _ in range(l):
        v = indel(v, 0)
    for _ in range(m):
        w = indel(w, 0)
    for _ in range(n):
        u = indel(u, 0)

    return align_score, v, w, u


if __name__ == "__main__":
    from dataimport import Score
    str1 = 'TWLNSACYGVNFRRLNPMNKTKWDCWTWVPMVMAAQYLCRIFIPVMDHWEFFGDWGLETWRLGIHDHVKIPNFRWSCELHIREHGHHFKTRFLKHNQFTQCYGLMPDPQFHRSYDVACQWEVTMSQGLMRFHRQNQIEKQRDRTSTYCMMTIGPGFTSNGYDPFVTITITPVQEPVENWFTPGGSMGFMIISRYMQMFFYLTRFSDMTYLVGVHCENYVCWNNVAKFLNGNLQGIFDQGERAYHQFVTWHSYSQYSRCSVGRYACEQAMSRVNSKMTWHWPIRDQGHEHFSEQYLSEKRNPPCNPRIGNAGQHFYEIHRIAHRVAMCNWAPQGQHPGGPTPHDVETCLWLWSLCLKGSDRGYVDRPWMFLADQLGEANLTLITMFHGCTRGCLMWFMDWEECVCSYSVVNPRCHGSEQWSVQNLGWRTCDTLISLWEPECDKHNTPPCLHWEFEDHPSQLRPVMMCDKYVQSIPTDAKWAWTYSKDFVISHWLIWTPIKLEECVFPQINRLWGTACNQGSQKIVIQNVWLRPSSFFQERSKCSDSSCILNVGGSNVNITGKETRTHVPILHMHEIDLISTASSGMRHNLILPHGMLMLHMNWHHSTRAMNPYSSLKLIPWTFQVCETDDRDQNVATHVADPCHKGEDQEIRCCKGGVDHQWKGDRMWMMCMPDMNYVKQDQAPSGTCEGACENYPADKDKCYMIFTIVFDYRRCTKKVCIWISGFPVDAFNLISIANAGFFCCWLEPTELKWRRTFYLGKGTQGWMCTFPHRNIIPVIICAGFGRWVQGEVPFRPVAQISAHSSDRRQGHHPPGTNMCHDYGDQYPIKRVGMQVEEDDGASYCDCAADWKLADMYEADHLSIGVIDFTDWIYPKNGGIWSEIIKSHFHWYHWETPQNTVGAFNTIVGINGSDMCIYHGNTQWEFGWCWKWLNHGHMRNQGPCHLGILEGRISKFAQVTSWWWQTKHDKDWSIEPYGRHWGEAGRPYTYNYCWMRWAIVYNHGNVISVELVPFMDEYPGKCNKEDVQFELFSPMQA'
    str2 = 'LWFKFLQCIFQYFKDQQETNCIWTFSPFSEHICQRVCQVYWNWNTPSSRTSDPRELFANSTIHNNRCGEWRYMFYHTRTLVQTAPLMKETLHSDGKHSMYCEQRHFFRSSYLIKVNYDVSHYLELYTFSEIPWKLTTHGWDGFSWFLLVNSCCTFDIDGKCGILSQCGMSRAFRTRQEDAYHFQTSLMHLHLHLHVQEGKHEKADLFAQFYNMLPMHGGTCGRNTEPSDLFDSATMNKYMAEHPASCKACPNVSKECFVYWWSHDFTKKHKLIEFSCGRDTGQTTQRTWNVDENEGGKWIWRFHYFMRAKALQIDPKFKPYWNEPRAIMRPGHVTAAPCICAQHSQNETAVCNRDQMHIHAIEFQQYHSRAFGEVQTWCDIGKENENDFIYEQHWWLVGGTEGMAGVIWKFVCARCRTQDCDFWKTCLTYSAQPMMKVYDTIFYVNSINPWEFEDHPSQCDKCVQSIPTDAKYAICGKFVISHWLYWTPQKFEECVHNNVRCAPMGNRLWGTACMVIQNVWLRPSMGSHFSCILNVGGSNINIQGKETWTHVPILHMHEIDLISTASSGMETCKPCFLSGPTIHMGFSYEIRAQPYSRDYFCMDWMQEADEVDHNRCETVQPTLPLLQQFEWKTSCMGQRWITIFCDHCQIVCFSTFFCVMPTFLPNTSILDKFYCIYLSISWTHYCNVHALGFIMRLHYSYMGWKEHKRMHAWDIGLDELWAQEGIQRAQLWCGDEFEVAKYPEWITEARTAIATRPWFHNCYIKPWWIREKHLWFGKESKLDHGHRGAMFTPVANDNTEWMHHWYMFCWAGSKNRLKRQIKEKLIFIIKFMITEFGLFLMIDYTQCYIAWMWAYTGIACYIDWEKCLKHDLTTTDLGCCVYRLFKWYEVRHRAPPQVNTRLPWSQIPMVAIQCNIVDECKEQWHFSYKASFVVEYLCPGCCTNGNRWQWYQVKETPFMYAFAASIFGFHHENLVVFITGSVTIPNGLFGCIAWTSPKPVQKTPASANTIIAYDKCILMG'
    penalty = 5
    score_matrix = Score().BLOSUM62
    m, back = scoring(penalty, score_matrix, str2, str1)
    with open('output.txt', 'w') as f:
        for row in m:
            f.write('\t'.join([str(x) for x in row]) + '\n')
