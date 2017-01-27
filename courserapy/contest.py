def chemical_reactions():
    f = open('input.txt', 'r')

    # overwrote builtin to emulate contest stdin
    def input():
        s = f.readline().strip('\n')
        if s == '':
            raise EOFError
        return s

    chemicals = set(map(int, input().split(' ')))
    reactions = dict()

    parse_rxn = lambda rxn: map(int, rxn.split('+'))

    while True:
        try:
            reactions[input()] = ''
        except EOFError:
            f.close()
            break

    for rxn in list(reactions):
        _ = list(filter(lambda s: s != '', rxn.split('->')))
        t = tuple(parse_rxn(_[0]))
        if t in reactions:
            [reactions[t].add(item) for item in parse_rxn(_[1])]
        else:
            reactions[t] = set(parse_rxn(_[1]))

        reactions.pop(rxn)

    while len(reactions) > 0:
        flag = False
        for p in list(reactions):
            if all([_ in chemicals for _ in p]):
                for _ in reactions[p]:
                    chemicals.add(_)
                reactions.pop(p)
                flag = True
                break

        if not flag:
            reactions.popitem()

    print(' '.join([str(x) for x in chemicals]))


def rna_folding():
    def foldr(str1, str2):

        return


def intron_detection(stream):
    # import data
    with open(stream, mode='r') as f:
        s, n, *reads = list(line.strip() for line in f.readlines())

    # score reads as introns or exons
        # ideas for scoring...
        #   % of read in guessed sequence
        #   index of read in guessed sequence

    # assemble some best guesses from string for actual DNA sequence
    # -> use best guesses of exons

    # refix scores of introns and exons

    # repeat until convergence

    # NOTES
        # initial guess for scoring will be difficult
        # maybe need to find most common substring among introns and exons
        # first?

        # or alternative method could be to create a splicing algo to knit
        # together potential reads with mismatches and test???
