import unittest

import __init__
import constants as cn
from depvar_group import IsolateMutationStatistics, DepvarGroup
import helpers
from study_context import StudyContext
from model_data_provider import ModelDataProvider
import util

import numpy as np
import pandas as pd


IGNORE_TEST = False
IS_PLOT = False


def dotestIsolateMutationStatistics(statistics):
  result = isinstance(statistics.ser_count, pd.Series)  \
      and (len(statistics.ser_count) > 5)
  return result


class TestIsolateMutationStatistics(unittest.TestCase):

  def setUp(self):
    constraints = [lambda r: r[cn.LINE] == cn.LINE_HA2]
    context = StudyContext(depvar=cn.RATE, 
        mutation_column=cn.GGENE_ID)
    self.provider = ModelDataProvider(context, 
        constraints=constraints)
    self.provider.do()
    self.statistics = IsolateMutationStatistics(self.provider)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertIsNone(self.statistics.ser_count)

  def testDo(self):
    if IGNORE_TEST:
      return
    self.statistics.do()
    self.assertTrue(dotestIsolateMutationStatistics(self.statistics))


class TestDepvarGroup(unittest.TestCase):

  def setUp(self):
    constraints = [lambda r: r[cn.LINE] == cn.LINE_HA2]
    self.group = DepvarGroup(cn.RATE, cn.GGENE_ID,
        constraints=constraints, is_plot=IS_PLOT)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(dotestIsolateMutationStatistics(
        self.group._statistics))

  def testMakeStatisticsDF(self):
    def test(lower, upper, steps):
      df = self.group.makeStatisticsDF(lower, upper, steps)
      self.assertTrue(helpers.isValidDataFrame(df,
          [cn.LOWER, cn.UPPER, cn.CENTER, cn.ISOLATES,
          cn.MUTATIONS, cn.FRACTION]))
      return df
    #
    for length in range(1, 5):
      df = test(-1, 1, length)
      self.assertEqual(len(df), length)

  def testMakeHeatmap(self):
    if IGNORE_TEST:
      return
    constraints = [lambda r: r[cn.LINE] == cn.LINE_HR2]
    for depvar in cn.DEPVARS:
      group = DepvarGroup(depvar, cn.GGENE_ID,
          constraints=constraints, is_plot=IS_PLOT)
      group.makeHeatmap(150)
    

if __name__ == '__main__':
    unittest.main()
