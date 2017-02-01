


if __name__ == "__main__":
    with open('input.txt', 'r') as f:
        genome1, genome2 = [line.strip() for line in f.readlines()]
    print(two_break(genome1, genome2))
