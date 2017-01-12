def greedy_sorting(permutation):
    # parse permutation into list of lists
    sign = {'+': '-', '-': '+'}
    permutation = permutation.strip('()').split(' ')
    P = []
    temp = []
    for item in permutation:
        temp.append(list(item))
    for i, item in enumerate(temp):
        j = i
        while P[0][j] !=

if __name__ == "__main__":
    greedy_sorting("(-3 +4 +1 +5 -2)")
