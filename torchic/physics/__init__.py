from torchic.physics.calibration import (
    py_BetheBloch,
    cluster_size_parametrisation,
)

from torchic.physics import ITS
from torchic.physics import simulations

import os

def try_import_root():
    try:
        import ROOT
        from ROOT import gInterpreter

        CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

        # Include headers and implementation
        gInterpreter.ProcessLine(f'#include "{CURRENT_DIR}/BetheBloch.hh"')

        # Import the class to make it accessible
        from ROOT import BetheBloch
        return BetheBloch

    except ImportError:
        print("ROOT not found. Functions will not be available.")
    except Exception as e:
        print(f"ROOT is available, but functions failed to compile: {e}")

    return None

BetheBloch = try_import_root()

__all__ = [
    'BetheBloch',
    'py_BetheBloch',
    'cluster_size_parametrisation',
    'ITS',
    'simulations',
]