'''
Written in Py35
Clever utilities, borrowed or self-wrote
'''


class Utils(object):

    @classmethod
    def test(cls, func, args=None, kwargs=None):
        '''
        canned test for profiling code
        iz naizz 4 large input tests
        '''
        import cProfile
        import pstats
        import io

        pr = cProfile.Profile()
        pr.enable()
        # no idea what this does...maybe CPU time in compiled C code...
        if args is not None and kwargs is not None:
            pr.runcall(func, *args, **kwargs)
        elif args is not None:
            pr.runcall(func, *args)
        else:
            pr.run(func)
        pr.disable()

        s = io.StringIO()
        ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
        ps.print_stats()
        print(s.getvalue())

    @classmethod
    def instance_size(cls, obj, seen=None):
        '''
        Borrowed from user:wissam on StackOverflow
        Returns the size of a Python object in memory

        nice b/c getsizeof doesn't account for size of underlying items
        '''
        from sys import getsizeof
        size = getsizeof(obj)
        if seen is None:
            seen = set()
        obj_hash = id(obj)
        if obj_hash in seen:
            return 0
        seen.add(obj_hash)

        # recurse
        if isinstance(obj, dict):
            size += sum([Utils.instance_size(_, seen)
                         for _ in list(obj.keys()) + list(obj.values())])
        elif hasattr(obj, '__dict__'):
            size += Utils.instance_size(obj.__dict__, seen)
        elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
            size += sum([Utils.instance_size(i, seen) for i in obj])

        # end recursion
        return size

    @classmethod
    def nested_dict(cls):
        '''
        Borrowed from user:andrew on StackOverflow
        Creates a defaultdict of arbitrary depth

        nice b/c it's a functional way of creating arbitrary depth
        '''
        from collections import defaultdict
        return defaultdict(Utils.nested_dict)


class Tree(object):
    '''
    very slow for collections above ~1000 items
    designed for ~100 items
    '''
    from sys import setrecursionlimit
    parent = None
    root = None
    left = None
    right = None
    setrecursionlimit(1000)

    def __init__(self, elem=None):
        if elem is not None and hasattr(elem, '__iter__') and not isinstance(elem, tuple):
            self.left = Tree()
            self.right = Tree()
            for e in elem:
                self.insert(e)
        else:
            if elem is not None:
                self.root = elem
                self.left = Tree()
                self.right = Tree()

    def __contains__(self, val):
        return bool(self[val])

    def __getitem__(self, key):
        if key == self.root:
            return self.root
        elif not self.left.isEmpty() and key < self.root:
            return self.left.__getitem__(key)
        elif not self.right.isEmpty() and key > self.root:
            return self.right.__getitem__(key)
        return

    def __setitem__(self, key, val):
        if key in self:
            self.remove(key)
            self.extend([key, val])

    def __iter__(self):
        if self.isEmpty():
            return
        if not self.left.isEmpty():
            for elem in self.left:
                yield elem
        yield self.root
        if not self.right.isEmpty():
            for elem in self.right:
                yield elem

    def hasParent(self):
        return bool(self.parent is not None)

    def isEmpty(self):
        return bool(self.root is None)

    def remove(self, key):
        try:
            self[key] = Tree()
        except KeyError as e:
            raise e

    def pop(self, key=None):
        if key is None:
            temp = self.root
            self.remove(self.root)
            return self.root
        else:
            temp = self.root
            self.remove(key)
            return temp

    def insert(self, leaf):
        if not self.isEmpty():
            if leaf < self.root:
                self.left.parent = self.root
                self.left.insert(leaf)
            elif leaf > self.root:
                self.right.parent = self.root
                self.right.insert(leaf)
            return
        else:
            self.root = leaf
            self.left = Tree()
            self.right = Tree()

    def extend(self, iterable):
        for elem in iterable:
            self.insert(elem)

    def union(self, tree):
        self.extend(item for item in tree)


if __name__ == '__main__':

    t = Tree(range(1000))
    for item in t:
        print(item)
