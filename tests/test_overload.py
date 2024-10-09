from torchic.utils.overload import overload

def test_function_overload():
    @overload
    def f(x: int) -> int:
        return x + 1

    @overload
    def f(x: int, y: int) -> int:
        return x + y

    print('f(1) = ',f(1))
    print('f(1, 2) = ',f(1, 2))
    assert f(1) == 2
    assert f(1, 2) == 3

if __name__ == '__main__':
    test_function_overload()
    print('Done')