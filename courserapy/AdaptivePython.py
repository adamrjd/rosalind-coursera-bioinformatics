
'''
def main():
    # generates a pretty star (': looks like this...
    #* . . . . . . * . . . . . . *
    #. * . . . . . * . . . . . * .
    #. . * . . . . * . . . . * . .
    #. . . * . . . * . . . * . . .
    #. . . . * . . * . . * . . . .
    #. . . . . * . * . * . . . . .
    #. . . . . . * * * . . . . . .
    #* * * * * * * * * * * * * * *
    #. . . . . . * * * . . . . . .
    #. . . . . * . * . * . . . . .
    #. . . . * . . * . . * . . . .
    #. . . * . . . * . . . * . . .
    #. . * . . . . * . . . . * . .
    #. * . . . . . * . . . . . * .
    #* . . . . . . * . . . . . . *
    n = int(input())
    matrix = list()
    is_star = lambda x, y: True if ((x == y)
                                    or (x == (n - 1) // 2)
                                    or ((x, y) == (n - 1 - y, n - 1 - x))
                                    or (y == (n - 1) // 2)) else False

    for _ in range(n):
        matrix.append(list('*' if is_star(_, __) else '.'
                        for __ in range(n)))

    print('\n'.join([' '.join(row) for row in matrix]))

'''

'''

class Matrix(object):
'''
# Interface to consume matrix data. If data is not a premade 2D array,
# the typeclass of string data to be parsed can be specified on
# initialization; by default it is int. Also, the print spacer is by
# default a tab.
'''
    matrix = None
    spacer = None
    N = None
    M = None

    def __init__(self, datum, spacer='\t', t=int):
        try:
            self.spacer = spacer

            assert isinstance(datum, list)
            self.matrix = datum
            self.M = len(datum)
            self.N = len(datum[0])
        except AssertionError:
            self.matrix = self._parse(t, datum)

    def __str__(self):
        return ''.join([_ + '\n' for _ in [self.spacer.join([str(_) for _ in row])
                                        for row in self.matrix]])

    def _parse(self, t, row_data):
        __ = [[t(x) for x in row.split(' ')]
            for row in filter(lambda _: _ != '', row_data.split('\n'))]
        self.M = len(__)
        self.N = len(__[0])
        return __

    def transpose(self):
        new = list(list() for _ in range(self.N))
        for _ in range(self.N):
            for __ in range(self.M):
                new[_].append(self.matrix[__][_])
        return new
        '''


'''class StringSet(object):
    size = None
    capacity = None
    strset = None
    m = None

    def __init__(self, args=None, s=8, m=0):
        self.size = 0
        self.m = m
        self.capacity = s
        self.strset = list(None for _ in range(self.capacity))
        if args is not None:
            for k in args:
                self.__setitem__(k, k)

    def __iter__(self):
        for __ in self.strset:
            if __ is not None:
                yield __

    def __getitem__(self, key):
        return self.strset[self.strhash(key)]

    def __setitem__(self, key, value):
        if (self.size + 1) / self.capacity >= 0.75:
            self.m += 1
            self.capacity = 8 ** (2 * self.m)
            self._rehash()
        _ = self.strhash(key)
        __ = self.strset[_]
        if __ is not None:
            if __ == key:
                # if adding something already in set
                if value is not None:
                    return
                # if removing key, i.e. value == None
                else:
                    self.strset[_] = None
                    return
            i = 1
            while __ is not None:
                _ = self.strhash(key, i)
                __ = self.strset[_]
                if __ == key:
                    return
                i += 1
            self.size += 1
            self.strset[_] = value
        else:
            self.size += 1
            self.strset[_] = value
        return

    def __contains__(self, key):
        return bool(self[key] is not None)

    def _rehash(self):
        __ = list(_ for _ in self)
        self.strset = [None for _ in range(self.capacity)]
        for key in __:
            self.__setitem__(key, key)

    def strhash(self, key, i=0):
        return (hash(key) + i ** 2) % self.capacity

    def add(self, key):
        self.__setitem__(key, key)
        return

    def remove(self, key):
        self.__setitem__(key, None)
        return


def main():
    f = open('input.txt', 'r')
    __ = 'a aa abC aa ac abc bcd a'
    text = list(map(lambda s: s.lower(), __.split(' ')))
    string_set = dict()
    for word in text:
        if word in string_set:
            string_set[word] += 1
        else:
            string_set[word] = 1
    for k in sorted(string_set, key=lambda k: string_set[k]):
        print(k + ' ' + str(string_set[k]))

main()'''
