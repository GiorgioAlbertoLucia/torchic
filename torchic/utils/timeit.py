'''
    Function decorator to time a function
'''

import time
from torchic.utils import TerminalColors as tc

def timeit(func):
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        print(tc.GREEN+'[INFO]: '+tc.RESET+f'{func.__name__} took {(end - start):.2f} seconds')
        return result
    return wrapper