import unittest
import random
import pandas as pd
from torchic.core.histogram import AxisSpec, build_TH1, build_TH2, fill_TH1, fill_TH2, build_efficiency, normalize_hist
from ROOT import TH1F, TH2F, TFile

class TestBuildHist(unittest.TestCase):

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)

    def setUp(self):
        self.data = pd.Series([1, 2, 3, 4, 5])
        self.data_2d = pd.DataFrame({'column_x': [1, 2, 3, 4, 5], 'column_y': [1, 2, 3, 4, 5]})
        self.axis_spec_x = AxisSpec(1, 5, 5, 'column_x', ';column_x;')
        self.axis_spec_y = AxisSpec(1, 5, 5, 'column_y', ';column_y;')

    def test_build_hist(self):
        hist = build_TH1(self.data, self.axis_spec_x)
        self.assertIsInstance(hist, TH1F)
        self.assertEqual(hist.GetEntries(), 5)
        self.assertEqual(hist.GetNbinsX(), self.axis_spec_x.nbins)
        self.assertEqual(hist.GetXaxis().GetXmin(), self.axis_spec_x.xmin)
        self.assertEqual(hist.GetXaxis().GetXmax(), self.axis_spec_x.xmax)

    def test_build_hist_2d(self):
        hist = build_TH2(self.data_2d['column_x'], self.data_2d['column_y'], self.axis_spec_x, self.axis_spec_y)
        self.assertIsInstance(hist, TH2F)
        self.assertEqual(hist.GetEntries(), 5)
        self.assertEqual(hist.GetNbinsX(), self.axis_spec_x.nbins)
        self.assertEqual(hist.GetNbinsY(), self.axis_spec_y.nbins)
        self.assertEqual(hist.GetXaxis().GetXmin(), self.axis_spec_x.xmin)
        self.assertEqual(hist.GetXaxis().GetXmax(), self.axis_spec_x.xmax)
        self.assertEqual(hist.GetYaxis().GetXmin(), self.axis_spec_y.xmin)
        self.assertEqual(hist.GetYaxis().GetXmax(), self.axis_spec_y.xmax)

    def test_fill_hist(self):
        hist = TH1F('hist', 'hist', 5, 1, 5)
        fill_TH1(self.data, hist)
        self.assertEqual(hist.GetEntries(), 5)

    def test_fill_hist_2d(self):
        hist = TH2F('hist', 'hist', 5, 1, 5, 5, 1, 5)
        fill_TH2(self.data_2d['column_x'], self.data_2d['column_y'], hist)
        self.assertEqual(hist.GetEntries(), 5)

    def test_build_efficiency(self):
        
        data_tot = [random.uniform(-0.5, 4.5) for _ in range(100)]
        data_sel = data_tot[:50]
        hist_tot = build_TH1(pd.Series(data_tot), AxisSpec(5, -0.5, 4.5, 'hist_tot', ';column;'))
        hist_sel = build_TH1(pd.Series(data_sel), AxisSpec(5, -0.5, 4.5, 'hist_sel', ';column;'))
        hist_eff = build_efficiency(hist_tot, hist_sel)
        output_file = TFile('data/test_build_efficiency.root', 'recreate')
        hist_tot.Write()
        hist_sel.Write()
        hist_eff.Write()
        output_file.Close()

    def test_normalize_hist(self):
        output_file = TFile('data/test_normalize_hist.root', 'recreate')
        data = [random.uniform(-0.5, 4.5) for _ in range(100)]
        hist = build_TH1(pd.Series(data), AxisSpec(5, -0.5, 4.5, 'hist', ';column;'))
        self.assertEqual(hist.Integral(), 100.)
        hist.Write()
        normalize_hist(hist)
        self.assertAlmostEqual(hist.Integral(), 1.)
        hist.Write('normalized')
        output_file.Close()
        
if __name__ == '__main__':
    unittest.main()