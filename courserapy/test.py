from operator import itemgetter
test_matrix = [list(_ for _ in range(__, __ * 2)) for __ in range(1, 6)]
print(test_matrix)
print(list(map(itemgetter(0), test_matrix)))
