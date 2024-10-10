import unittest
import pandas as pd
from torchic.core.histogram import AxisSpec
from torchic.core.dataset import Dataset

class TestDataset(unittest.TestCase):

    def setUp(self):
        self.data = pd.DataFrame({
            'column_x': [1, 2, 3, 4, 5],
            'column_y': [5, 4, 3, 2, 1]
        })
        self.dataset = Dataset(self.data)

    def test_data_property(self):
        self.assertTrue(self.dataset.data.equals(self.data))

    def test_add_subset(self):
        self.dataset.add_subset('subset', self.data['column_x'] > 2)
        self.assertEqual(len(self.dataset.subsets['subset']), 3)

    def test_build_hist(self):
        axis_spec_x = AxisSpec(1, 5, 5, 'column_x', ';column_x;')
        hist = self.dataset.build_hist('column_x', axis_spec_x)
        self.assertEqual(hist.GetEntries(), 5)

    def test_build_hist_2d(self):
        axis_spec_x = AxisSpec(1, 5, 5, 'column_x', ';column_x;column_y;')
        axis_spec_y = AxisSpec(1, 5, 5, 'column_y', ';column_x;column_y;')
        hist = self.dataset.build_hist('column_x', 'column_y', axis_spec_x, axis_spec_y)
        self.assertEqual(hist.GetEntries(), 5)

    def test_build_hist_subset(self):
        self.dataset.add_subset('subset', self.data['column_x'] > 2)
        axis_spec_x = AxisSpec(1, 5, 5, 'column_x', ';column_x;')
        hist = self.dataset['subset'].build_hist('column_x', axis_spec_x)
        self.assertEqual(hist.GetEntries(), 3)
    
    def test_build_hist_subset_alternative(self):
        self.dataset.add_subset('subset', self.data['column_x'] > 2)
        axis_spec_x = AxisSpec(1, 5, 5, 'column_x', ';column_x;')
        hist = self.dataset.build_hist('column_x', axis_spec_x, subset='subset')
        self.assertEqual(hist.GetEntries(), 3)

    def test_query(self):
        self.dataset.query('column_x > 2')
        self.assertEqual(len(self.dataset.data), 3)

if __name__ == '__main__':
    unittest.main()