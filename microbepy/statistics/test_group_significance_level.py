"""Tests for simple_statistics"""

import microbepy_init
import util
from group_significance_level import GroupSignificanceLevel
from model_data_provider import ModelDataProvider
from mutation_context import MutationContext
import util_data as ud
import constants as cn
import helpers

import numpy as np
import os
import pandas as pd
import unittest

IGNORE_TEST = False

MUTATIONS = ['DVU1295', 'MMP1314']

########################################
class TestGroupSignificanceLevel(unittest.TestCase):

  def setUp(self):
    self.cls = GroupSignificanceLevel
    self.provider = ModelDataProvider(
        MutationContext(cn.YIELD, cn.GENE_ID))
    self.provider.do()
    self.significance = GroupSignificanceLevel(self.provider,
        MUTATIONS)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(helpers.isValidDataFrame(
        self.significance.df_group,
        self.significance.df_group.columns))
    for sl in [self.significance.sl_min, self.significance.sl_max]:
      self.assertEqual(self.significance._mutations,
          sl.mutations)
    self.assertLess(self.significance.sl_min.avg,
        self.significance.sl_max.avg)
    self.assertFalse(self.significance.sl_min.group ==
        self.significance.sl_max.group)

  def testMakeGroupStatisticDF(self):
    if IGNORE_TEST:
      return
    df = self.significance.makeGroupStatisticDF()
    self.assertTrue(helpers.isValidDataFrame(df,
        [cn.GROUP, cn.AVG, cn.STD, cn.COUNT]))

  def testCalcMinMaxSLTStatBasic(self):
    if IGNORE_TEST:
      return
    self.significance.calcMinMaxSLTStat()
    sl_min = self.significance.sl_min
    sl_max = self.significance.sl_max
    trues = [(v <= 1.0) and (v >= 0) 
          for v in [sl_min.sl_tstat, sl_max.sl_tstat]]
    self.assertTrue(all(trues))

  def testCalcMinMaxSLTStatData(self):
    if IGNORE_TEST:
      return
    """
    Test the calculations. Do on
    random mutations and then on a mutation combination that
    looks like it should be significant.
    """
    mutations = ['DVU1295', 'MMP1314']
    significance = self.cls(self.provider, mutations)
    significance.calcMinMaxSLTStat()
    sl_min = significance.sl_min
    sl_max = significance.sl_max
    self.assertLess(sl_min.sl_tstat, 0.5)
    self.assertLess(sl_max.sl_tstat, 0.02)

  def testCalcMinMaxSLResample(self):
    if IGNORE_TEST:
      return
    """
    Test the calculations. Do on
    random mutations and then on a mutation combination that
    looks like it should be significant.
    """
    def test(sl):
      self.assertGreaterEqual(sl, 0)
      self.assertLessEqual(sl, 1.0)
    MUTATION_SETS = [
        ['DVU1295', 'MMP1314'],
        ['DVU2451', 'MMP1255']
        ]
    for mutations in MUTATION_SETS:
      significance = GroupSignificanceLevel(self.provider,
          mutations)
      significance.calcMinMaxSLResample(num_replications=int(1e4))
      test(significance.sl_min.sl_resample)
      test(significance.sl_max.sl_resample)

  def testGetStatistics(self):
    self.significance.calcMinMaxSLTStat()
    self.significance.calcMinMaxSLResample()
    sl_min, sl_max = self.significance.getStatistics()
    self.assertFalse(np.isnan(sl_min))
    self.assertFalse(np.isnan(sl_max))
    
    
if __name__ == '__main__':
    unittest.main()
