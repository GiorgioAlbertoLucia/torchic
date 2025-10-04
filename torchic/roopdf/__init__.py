import os

def try_import_root():
    try:
        import ROOT
        from ROOT import gInterpreter

        CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
        pdf_dir = os.path.join(CURRENT_DIR, 'RooCustomPdfs')

        # Include headers and implementation
        gInterpreter.ProcessLine(f'#include "{pdf_dir}/RooGausExp.hh"')
        gInterpreter.ProcessLine(f'#include "{pdf_dir}/RooSillPdf.hh"')
        gInterpreter.ProcessLine(f'#include "{pdf_dir}/RooGausDExp.hh"')

        # Import the class to make it accessible
        from ROOT import RooGausExp, RooSillPdf, RooGausDExp
        return RooGausExp, RooSillPdf, RooGausDExp

    except ImportError:
        print("ROOT not found. Functions will not be available.")
    except Exception as e:
        print(f"ROOT is available, but functions failed to compile: {e}")

    return None, None

RooGausExp, RooSillPdf, RooGausDExp = try_import_root()

__all__ = [
    'RooGausExp',
    'RooSillPdf',
    'RooGausDExp',
]