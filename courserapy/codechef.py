
'''
# Ordering the Soldiers
# ORDERS

f = open('input.txt', 'r')


def input():
    __ = f.readline().strip('\n')
    return __ if __ != '' else None


def main():
    def run():
        n = int(input())
        order = list(map(int, input().split(' ')))
        soldiers = list(str(i) for i in range(1, n + 1))
        for i, o in enumerate(order[::-1]):
            if o > 0:
                soldiers.insert(n - i - 1, soldiers.pop(n - i - 1 - o))
        print(' '.join(soldiers))
    num_cases = int(input())
    for __ in range(num_cases):
        run()

main()
'''
'''
#Nothing in Common
# NOTINCOM


def main():
    N = int(input())
    for _ in range(N):
        run()


def run():
    __ = '3 4\n1 2 3\n3 4 5 6'.split('\n')
    parse_and_map = lambda: list(map(int, __[0].split(' ')))
    (N, M) = parse_and_map()
    __.pop(0)
    A, B = dict(), dict()
    for _, d in enumerate([A, B]):
        for k in parse_and_map():
            if k in d:
                d[k] += 1
            else:
                d[k] = 1
        __.pop(0)
    print(sum(list(min(A[k], B[k]) for k in A if k in B)))
run()
'''
'''
Extremely Large Inputs

def check(num, k):
    return 1 if num % k == 0 else 0

n, k = map(int, input().split(' '))

results = []
for i in range(n):
    results.append(check(int(input()), k))

print(sum(results))'''

'''
Digit Longest Increasing Subsequences 2
LISDIGIT



def solve(line):
    answer = []

    return int(''.join([str(x) for x in answer]))
n = int(input())
results = dict()  # more space efficient than list
for i in range(n):
    input()  # first line is number of digits...can be inferred
    results[i] = solve(list(map(int, input().split(' '))))

for i in range(n):
    print(results[i])


print(bin(5))
print(bin(3))

print(bin(10))

print(bin(196613))
print(bin(655370))


print(bin(9999 | (9999 << (1 << (1 << (1 << 1))))))

print(bin(9999))

'''
