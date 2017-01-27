n, m = list(map(int, input()))
elem_with_spaces = lambda e: ' ' * (4 - len(str(e))) + str(e)
matrix = []

for row in range(n):
    r = range(m * row, m * (row + 1))
    if row % 2 == 1:
        r = reversed(r)
    matrix.append(''.join([elem_with_spaces(e) for e in r]))

print(''.join([elem + '\n' for elem in matrix]))
