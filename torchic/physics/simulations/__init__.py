import os
def try_import_root():
    try:
        import ROOT
        from ROOT import gInterpreter

        CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

        # Include headers and implementation
        gInterpreter.ProcessLine(f'#include "{CURRENT_DIR}/ExponentialDecaySimulation.cxx"')
        gInterpreter.ProcessLine(f'#include "{CURRENT_DIR}/TwoBodyDecaySimulation.cxx"')
        gInterpreter.ProcessLine(f'#include "{CURRENT_DIR}/test.cxx"')

        # Import the class to make it accessible
        from ROOT import RunExponentialDecaySimulation, RunTwoBodyDecaySimulation, hello
        return RunExponentialDecaySimulation, RunTwoBodyDecaySimulation, hello

    except ImportError:
        print("ROOT not found. Functions will not be available.")
    except Exception as e:
        print(f"ROOT is available, but functions failed to compile: {e}")

    return None, None, None

RunExponentialDecaySimulation, RunTwoBodyDecaySimulation, hello = try_import_root()

__all__ = [
    'RunExponentialDecaySimulation',
    'RunTwoBodyDecaySimulation',
    'hello',
]