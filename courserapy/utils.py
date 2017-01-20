'''
Written in Py35
Clever utilities that I either borrowed or wrote myself

For now, it's just borrowed code ;____; 
'''


def nested_dict():
    '''
    Borrowed from user:andrew on StackOverflow
    Creates a defaultdict of arbitrary depth
    '''
    from collections import defaultdict
    return defaultdict(nested_dict)


def get_size(obj, seen=None):
    '''
    Borrowed from user:wissam on StackOverflow
    Returns the size of a Python object in memory
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
        size += sum([get_size(_, seen)
                     for _ in list(obj.keys()) + list(obj.values())])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])

    # end recursion
    return size

if __name__ == '__main__':
    from time import time

    t0 = time()
    test_dict = {(i, j): 1 for i, j in zip(range(1000), range(1000))}
    t_dict = time() - t0
    s_dict = get_size(test_dict)

    t0 = time()
    test_list = [[1 for _ in range(1000)] for __ in range(1000)]
    t_list = time() - t0
    s_list = get_size(test_list)

    print('Dictionary size and time complexity...')
    print(s_dict)
    print(t_dict)
    print()
    print('List space and time complexity...')
    print(s_list)
    print(t_list)
