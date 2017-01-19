def chromosome_to_cycle(chromosome):
    nodes = [0 for _ in range(len(chromosome) * 2)]
    for j, i in enumerate(chromosome):
        i = int(i)
        if i > 0:
            nodes[2 * j] = 2 * i - 1
            nodes[2 * j + 1] = 2 * i
        else:
            nodes[2 * j] = -2 * i
            nodes[2 * j + 1] = -2 * i - 1
    return nodes


def cycle_to_chromosome(nodes):
    chromosome = [0 for _ in range(len(nodes) // 2)]
    for j in range(len(nodes) // 2):
        if nodes[2 * j] < nodes[2 * j + 1]:
            chromosome[j] = nodes[2 * j + 1] // 2
        else:
            chromosome[j] = - nodes[2 * j] // 2
    return '(' + ' '.join(
        ['+' + str(x) if x > 0 else str(x)
         for x in chromosome]) + ')'


def colored_edges()

if __name__ == "__main__":
    pass
