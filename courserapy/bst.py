<<<<<<< HEAD
def main():
    n = int(input())
    matrix = [list(map(int, input().split(' '))) for _ in range(n)]
    flag = True

    for i, row in enumerate(matrix):
        if row != list(range(i, n + i)):
            flag = False
            break

    if flag:
        print('YES')
    else:
        print('NO')

print([1, 2, 3] == list(range(1, 4)))
=======
n, m = list(map(int, input()))
elem_with_spaces = lambda e: ' ' * (4 - len(str(e))) + str(e)
matrix = []

for row in range(n):
    r = range(m * row, m * (row + 1))
    if row % 2 == 1:
        r = reversed(r)
    matrix.append(''.join([elem_with_spaces(e) for e in r]))

print(''.join([elem + '\n' for elem in matrix]))
>>>>>>> 764fdc1829d88311c1202532396aea8f0513b53d
