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
'''

print(bin(5))
print(bin(3))

print(bin(10))

print(bin(196613))
print(bin(655370))


print(bin(9999 | (9999 << (1 << (1 << (1 << 1))))))

print(bin(9999))
