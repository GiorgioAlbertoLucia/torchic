import unittest
import pandas as pd
from torchic.core.histogram import AxisSpec, build_hist, fill_hist
from ROOT import TH1F, TH2F

class TestBuildHist(unittest.TestCase):

    def setUp(self):
        self.data = pd.Series([1, 2, 3, 4, 5])
        self.data_2d = pd.DataFrame({'column_x': [1, 2, 3, 4, 5], 'column_y': [1, 2, 3, 4, 5]})
        self.axis_spec_x = AxisSpec(1, 5, 5, 'column_x', ';column_x;')
        self.axis_spec_y = AxisSpec(1, 5, 5, 'column_y', ';column_y;')

    def test_build_hist(self):
        hist = build_hist(self.data, self.axis_spec_x)
        self.assertIsInstance(hist, TH1F)
        self.assertEqual(hist.GetEntries(), 5)
        self.assertEqual(hist.GetNbinsX(), self.axis_spec_x.nbins)
        self.assertEqual(hist.GetXaxis().GetXmin(), self.axis_spec_x.xmin)
        self.assertEqual(hist.GetXaxis().GetXmax(), self.axis_spec_x.xmax)

    def test_build_hist_2d(self):
        hist = build_hist(self.data_2d['column_x'], self.data_2d['column_y'], self.axis_spec_x, self.axis_spec_y)
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
        fill_hist(self.data, hist)
        self.assertEqual(hist.GetEntries(), 5)

    def test_fill_hist_2d(self):
        hist = TH2F('hist', 'hist', 5, 1, 5, 5, 1, 5)
        fill_hist(self.data_2d['column_x'], self.data_2d['column_y'], hist)
        self.assertEqual(hist.GetEntries(), 5)

if __name__ == '__main__':
    unittest.main()