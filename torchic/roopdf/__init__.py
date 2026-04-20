import os

class RooPdf:
    @staticmethod
    def try_import_roopdf():
        try:
            import ROOT
            from ROOT import gInterpreter

            CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
            pdf_dir = os.path.join(CURRENT_DIR, 'RooCustomPdfs')

            # Include headers and implementation
            gInterpreter.ProcessLine(f'#include "{pdf_dir}/RooGausExp.hh"')
            gInterpreter.ProcessLine(f'#include "{pdf_dir}/RooSillPdf.hh"')
            gInterpreter.ProcessLine(f'#include "{pdf_dir}/RooSillGeneralizedPdf.hh"')
            gInterpreter.ProcessLine(f'#include "{pdf_dir}/RooSillGeneralizedKstarPdf.hh"')
            gInterpreter.ProcessLine(f'#include "{pdf_dir}/RooGausDExp.hh"')

            # Import the class to make it accessible
            from ROOT import RooGausExp, RooSillPdf, RooSillGeneralizedPdf, RooSillGeneralizedKstarPdf, RooGausDExp
            return RooGausExp, RooSillPdf, RooSillGeneralizedPdf, RooSillGeneralizedKstarPdf, RooGausDExp

        except ImportError:
            print("ROOT not found. Functions will not be available.")
        except Exception as e:
            print(f"ROOT is available, but functions failed to compile: {e}")

        return None, None, None, None, None

RooGausExp, RooSillPdf, RooSillGeneralizedPdf, RooSillGeneralizedKstarPdf, RooGausDExp = RooPdf.try_import_roopdf()

__all__ = [
    'RooGausExp',
    'RooSillPdf',
    'RooSillGeneralizedPdf',
    'RooSillGeneralizedKstarPdf',
    'RooGausDExp',
]
