import os
from ROOT import gInterpreter
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

HELLO_PATH = os.path.join(CURRENT_DIR, 'test.cxx')
gInterpreter.ProcessLine(f'#include "{HELLO_PATH}"')
from ROOT import hello

EXPDECAYSIM_PATH = os.path.join(CURRENT_DIR, 'ExponentialDecaySimulation.cxx')
gInterpreter.ProcessLine(f'#include "{EXPDECAYSIM_PATH}"')
from ROOT import RunExponentialDecaySimulation

TOWBODYDECAYSIM_PATH = os.path.join(CURRENT_DIR, 'TwoBodyDecaySimulation.cxx')
gInterpreter.ProcessLine(f'#include "{TOWBODYDECAYSIM_PATH}"')
from ROOT import RunTwoBodyDecaySimulation