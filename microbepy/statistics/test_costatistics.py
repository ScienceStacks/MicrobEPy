"""Tests for simple_statistics"""

import microbepy_init
import util
from costatistics import CoStatistics
import costatistics as co
import constants as cn
import helpers
import util_data

import numpy as np
import os
import pandas as pd
import unittest

IGNORE_TEST = False
DF_DATA = util_data.makeIsolateData(is_separate_species=False)

########################################
class TestCoStatistics(unittest.TestCase):

  def setUp(self):
    if IGNORE_TEST:
      return
    self.statistics = CoStatistics(df_data=DF_DATA)

  def testMakeCountDF(self):
    if IGNORE_TEST:
      return
    def test(row_label, for_label, checkFunc, min_rows=1):
      """
      Tests if counts are collected for the keys.
      :param str row_label: how statistics are grouped
      :param str for_label: what is counted for each grouping
      :param function checkFunc: boolean function of the list of row_label values
      """
      df = self.statistics.makeCountDF(row_label=row_label, 
          for_label=for_label)
      trues = checkFunc([x for x in df[row_label].tolist()])
      self.assertTrue(trues)
      expected_columns = [row_label, cn.COUNT]
      self.assertTrue(helpers.isValidDataFrame(df, 
          expected_columns, min_rows=min_rows))

    checkFunc = lambda xs: any([x[0] == cn.SPECIES_MIX_DVH for x in xs])  \
        and any([x[0] == cn.SPECIES_MIX_DVH for x in xs])
    test(cn.KEY_MUTATION, cn.KEY_ISOLATE, checkFunc)
    test(cn.GENE_ID, cn.KEY_ISOLATE, checkFunc)
    test(cn.GENE_ID, cn.KEY_CULTURE, checkFunc)

  def testMakeParitionedCountDF2(self):
    df = self.statistics.makeCountDF(row_label=cn.KEY_MUTATION,
        for_label=cn.KEY_ISOLATE, col_label=cn.LINE)
    self.assertTrue(helpers.isValidDataFrame(df, df.columns))

    
if __name__ == '__main__':
    unittest.main()
