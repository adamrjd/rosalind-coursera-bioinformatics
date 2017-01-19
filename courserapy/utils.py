def nested_dict():
    '''
    Borrowed from user:andrew on StackOverflow
    Creates a defaultdict of arbitrary depth
    '''
    from collections import defaultdict
    return defaultdict(nested_dict)
