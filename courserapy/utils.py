'''
Written in Py35
Recipes I like or made
'''


class Utils(object):

    @staticmethod
    def groupN(iterable, n):
        '''s -> (s0,s1,s2,..., sn-1), (sn, sn+1, sn+2,...,s2n-1), ...'''
        return zip(*[iter(iterable)] * n)

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


if __name__ == '__main__':
    pass
