'''
    Function decorator to time a function
'''

import time

def timeit(func):
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        print(f'{func.__name__} took {(end - start):.2f} seconds')
        return result
    return wrapper