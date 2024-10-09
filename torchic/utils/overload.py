'''
    Classes and functions used for function overloading
'''
import inspect

class overload:
    # This class will act as the decorator directly
    functions = {}

    def __init__(self, func):
        self.func = func
        self.register(func)

    def register(self, func):
        # Get the function signature (number and types of arguments)
        signature = inspect.signature(func)
        param_types = tuple(param.annotation for param in signature.parameters.values())

        if func.__name__ not in overload.functions:
            overload.functions[func.__name__] = {}

        # Store function by its name and argument types
        overload.functions[func.__name__][param_types] = func

    def __call__(self, *args, **kwargs):
        func_key = tuple(type(arg) for arg in args)
        func_name = self.func.__name__

        if func_key in overload.functions[func_name]:
            # Dispatch to the appropriate overloaded function
            return overload.functions[func_name][func_key](*args, **kwargs)
        else:
            raise TypeError(f"No matching function for arguments {func_key}")