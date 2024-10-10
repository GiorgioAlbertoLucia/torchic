'''
    Freebooted from pythonlangutil
'''
def signature(*types):
    def func(f):
        def inner_func(callingObj, *args, **kwargs):
            if callingObj is None:
                return f(*args, **kwargs)
            else:
                return f(callingObj, *args, **kwargs)
        inner_func.signature = types
        return inner_func
    return func
        
class overload(object):
    def __init__(self, func):
        self.owner = None
        self.signatures = []
        self.methods = []
        self.methods.append(func)
        self.signatures.append(func.signature)
        
    def __get__(self, owner, ownerType=None):
        self.owner = owner or self
        return self
    
    def __call__(self, *args, **kwargs):
        '''
            To each version of the overloaded method, a signature is assigned. 
            The signature is a list of types of the accepted arguments.
            The __call__ method checks the signature of the arguments passed to the method and
            calls the correct version of the overloaded method. 
            NOTE: The signature does not consider kwargs.
        '''
        signature = []
        for arg in args:
            signature.append(arg.__class__.__name__)
        signature = tuple(signature)
        if signature in self.signatures:
            index = self.signatures.index(signature)
        else:
            raise Exception("There is no overload for this method with this signature.")
        return self.methods[index](self.owner, *args, **kwargs)
    
    def overload(self, func):
        self.methods.append(func)
        self.signatures.append(func.signature)
        return self