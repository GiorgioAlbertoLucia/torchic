import os
from ROOT import gInterpreter
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
EXPDECAY_DIR = os.path.join(CURRENT_DIR, 'ExponentialDecaySimulation.cpp')
gInterpreter.ProcessLine(f'#include "{EXPDECAY_DIR}"')
TWOBODYDECAY_DIR = os.path.join(CURRENT_DIR, 'TwoBodyDecaySimulation.cpp')
gInterpreter.ProcessLine(f'#include {TWOBODYDECAY_DIR}')
from ROOT import RunExponentialDecaySimulation, RunTwoBodyDecaySimulation