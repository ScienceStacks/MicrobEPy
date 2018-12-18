from dataframe_sorter import DataframeSorter

import unittest
import numpy as np
import pandas as pd


IGNORE_TEST = False
DF = pd.DataFrame({
    'b': [np.nan, 0, 0, 1],
    'd': [1, 1, 1, 1],
    'c': [1, 1, 1, 0],
    'a': [0, 0, 0, 0],
    }).transpose()
#
DF_SYMMETRIC = pd.DataFrame({
    'a': [-1, -1, -1, -1],
    'b': [-1, 0, 0, -1],
    'c': [-1, 0, 0, -1],
    'd': [-1, -1, -1, -1],
    })
DF_SYMMETRIC.index = DF_SYMMETRIC.columns
#
x = np.nan
DF_NANS = pd.DataFrame({
    'a': [1, x, x, 1],
    'b': [x, x, x, x],
    'c': [x, 1, x, x],
    'd': [x, x, x, x],
    })
DF_NANS_NONNAN = ['a', 'c']


class TestDataframeSorter(unittest.TestCase):

  def setUp(self):
    self.sorter = DataframeSorter(DF)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(isinstance(self.sorter._df, pd.DataFrame))

  def orderingTest(self, df_result, df_reference):
    self.assertEqual(len(df_reference), len(df_result))
    self.assertEqual(set(df_reference.columns), set(df_result))

  def testOrderRows(self):
    if IGNORE_TEST:
      return
    df_result = self.sorter.orderRows(max_clusters=2)
    self.orderingTest(df_result, DF)
    trues = [v == t for v, t in 
        zip(DF.index.values, ['a', 'b', 'c', 'd'])]
    self.assertTrue(all(trues))

  def testOrderRows(self):
    if IGNORE_TEST:
      return
    sorter = DataframeSorter(DF)
    df_result = sorter.orderColumns(is_convert_zero_to_nan=False)
    self.assertEqual(set(DF.columns), set(df_result.columns))
    self.assertEqual(set(DF.index), set(df_result.index))
    self.orderingTest(df_result.transpose(), DF.transpose())

  def testOrderBoth(self):
    if IGNORE_TEST:
      return
    sorter = DataframeSorter(DF_SYMMETRIC)
    df = sorter.orderBoth(is_convert_zero_to_nan=False,
        max_clusters=2)
    self.assertTrue(df.equals(df.transpose()))

  def testDeleteNanRowsAndColumns(self):
    if IGNORE_TEST:
      return
    sorter = DataframeSorter(DF_NANS)
    df = sorter.deleteNanRowsAndColumns()
    self.assertEqual(set(DF_NANS_NONNAN), set(DF_NANS.columns))
    self.assertGreater(len(DF_NANS), len(df))
    

if __name__ == '__main__':
    unittest.main()
