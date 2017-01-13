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


def two_break(p, q):
    from math import fabs
    cycles = lambda perm: len(perm.split(")("))
    breaks = lambda perm: breakpoint_count(perm.replace(")(", " "))
    return int(fabs((cycles(q) - breaks(q))))


def chromosome_to_cycle(chromosome):
    nodes = [0 for _ in range(2 * len(chromosome))]
    for j, i in enumerate(chromosome):
        if i > 0:
            nodes[2 * j - 1] = 2 * i - 1
            nodes[2 * j] = 2 * i
        else:
            nodes[2 * j - 1] = -2 * i
            nodes[2 * j] = -2 * i - 1
    return [nodes.pop()] + nodes


def cycle_to_chromosome(nodes):
    chromosomes = list()
    for j in range(len(nodes) // 2):
        if nodes[2 * j] < nodes[2 * j - 1]:
            chromosomes.append(nodes[2 * j - 1] // 2)
        else:
            chromosomes.append(- nodes[2 * j] // 2)
    return chromosomes

if __name__ == "__main__":
    print(6 % 5)
