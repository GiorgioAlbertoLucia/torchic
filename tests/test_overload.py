from torchic.utils.overload import overload, signature

def test_function_overload():
    @overload
    @signature('int')
    def f(x: int) -> int:
        return x + 1

    @f.overload
    @signature('int', 'int')
    def f(x: int, y: int) -> int:
        return x + y

    print('f(1) = ',f(1))
    print('f(1, 2) = ',f(1, 2))
    assert f(1) == 2
    assert f(1, 2) == 3

class MyClass:
    
    @overload
    @signature('int')
    def foo(self, a: int):
        return f"Method with integer called with {a}"

    @foo.overload
    @signature('str')
    def foo(self, a: str):
        return f"Method with string called with {a}"

    @foo.overload
    @signature('int', 'int')
    def foo(self, a: int, b: int):
        return f"Method with two integers called: {a}, {b}"
    
def test_class_overload():
    obj = MyClass()
    print('obj.foo(5) = ',obj.foo(5))
    assert obj.foo(5) == "Method with integer called with 5"
    print('obj.foo("hello") = ',obj.foo("hello"))
    assert obj.foo("hello") == "Method with string called with hello"
    print('obj.foo(5, 10) = ',obj.foo(5, 10))
    assert obj.foo(5, 10) == "Method with two integers called: 5, 10"

if __name__ == '__main__':
    test_function_overload()
    test_class_overload()
    print('Done')