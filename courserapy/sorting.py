def parse(permutation):
    return [int(_) for _ in permutation.strip('()').split(' ')]


def greedy_permutation_sort(p):
    permutation = parse(p)
    length = len(permutation)
    sign = {'+': '-', '-': '+'}
    reverse = lambda ls: list(reversed([(sign[l[0]], l[1]) for l in ls]))
    history = []
    if permutation[0][0] == '+' and permutation[0][1] == 1:
        history.append(permutation)
    c = 0
    while True:
        for i in range(length):
            p = permutation[i]
            if p[0] == '+' and p[1] == i + 1:
                c = i + 1
            else:
                j = i
                while p[1] != c + 1:
                    j += 1
                    p = permutation[j]
                permutation = permutation[
                    :i] + reverse(permutation[i:j + 1]) + permutation[j + 1:]
                history.append(permutation)
                break
        if c == length:
            break
    return history


def breakpoint_count(p):
    from math import fabs
    permutation = parse(p)
    permutation.insert(0, 0)
    permutation.append(len(permutation))
    points = 0
    for i in range(1, len(permutation)):
        perm = permutation[i - 1]
        next_perm = permutation[i]
        if fabs(fabs(perm) - fabs(next_perm)) != 1:
            points += 1
        else:
            if (perm > 0 and next_perm < 0) or (perm < 0 and next_perm > 0):
                points += 1
            elif (perm > 0 and next_perm > 0 and perm > next_perm):
                points += 1
            elif (perm < 0 and next_perm < 0 and fabs(perm) < fabs(next_perm)):
                points += 1
    return points


def chromosome_to_cycle(chromosome):
    try:
        assert isinstance(chromosome, list)
    except AssertionError:
        # parse chromosome string if necessary
        chromosome = [int(x) for x in chromosome.strip('()').split(' ')]

    cycle = [0 for _ in range(len(chromosome) * 2)]
    for j, i in enumerate(chromosome):
        i = int(i)
        if i > 0:
            cycle[2 * j] = 2 * i - 1
            cycle[2 * j + 1] = 2 * i
        else:
            cycle[2 * j] = -2 * i
            cycle[2 * j + 1] = -2 * i - 1
    return cycle


def cycle_to_chromosome(cycle):
    try:
        assert isinstance(cycle, list)
    except AssertionError:
        # parse cycle string if necessary
        cycle = [int(x) for x in cycle.strip('()').split(' ')]

    chromosome = [0 for _ in range(len(cycle) // 2)]
    for j in range(len(cycle) // 2):
        if cycle[2 * j] < cycle[2 * j + 1]:
            chromosome[j] = cycle[2 * j + 1] // 2
        else:
            chromosome[j] = - cycle[2 * j] // 2
    return '(' + ' '.join(
        ['+' + str(x) if x > 0 else str(x)
         for x in chromosome]) + ')'


def colored_edges(genome):
    try:
        assert isinstance(genome, list)
    except AssertionError:
        # parse genome string if necessary
        from re import split
        genome = [[int(x) for x in _.split(' ')]
                  for _ in filter(lambda _: _ != '', split(r'[)(]', genome.strip('()')))]

    edges = set()
    for chromosome in genome:
        cycle = chromosome_to_cycle(chromosome)
        for j in range(len(chromosome)):
            edges.add((cycle[2 * j - 1], cycle[2 * j]))
    return list(sorted(edges, key=lambda _: _[0]))


def genome_to_edges(genome):
    return colored_edges(genome)

'''def edges_to_genome(edges):
    try:
        assert isinstance(edges, list)
    except AssertionError:
        # parse edges string if necessary
        edges = [[int(x) for x in _.strip('()').split(', ')]
                 for _ in edges.split('), ')]

    genome = set()
    for cycle in edges:
        nodes = [_ for _ in range(cycle[0], cycle[1] + 1)]
        chromosome = cycle_to_chromosome(nodes)
        genome.add(chromosome)
    return genome'''


def two_break(p, q):
    blocks = lambda _: len([__ for __ in colored_edges(_) if __[1] < __[0]])
    return blocks(q) - breakpoint_count(p)


if __name__ == "__main__":
    with open('input.txt', 'r') as f:
        genome1, genome2 = [line.strip() for line in f.readlines()]
    print(two_break(genome1, genome2))
