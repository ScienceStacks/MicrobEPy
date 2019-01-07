import microbepy_init
import helpers
import isolate_regression as ir
import isolate_model as rm
import regression_plot as rp
from isolate import Isolate
import util

import copy
import collections
import constants as cn
import numpy as np
import pandas as pd
import os
import unittest


DIR_SEPARATOR = "/"
IGNORE_TEST = False
DEPVAR = cn.RATE
ISOLATES = ["HA2.152.09.01.D.CI", "HA2.152.09.01.M.CI"]
ISOLATES2 = ["HA2.152.09.02.D.CI", "HA2.152.09.02.M.CI"]
LINE = "HA2"
SIZE_OF_CROSS_RESULT = 4  # Number of elements in a cross validation

class TestFunctions(unittest.TestCase):

  def testMakeAvgAndStdDF(self):
    if IGNORE_TEST:
      return
    DATA = range(10)
    COL_A = 'col_a'
    COL_B = 'col_b'
    data = range(4)
    df = pd.DataFrame({
      COL_A: ['a', 'a', 'b', 'b'],
      COL_B: data,
      })
    df_avg_std = ir.makeAvgAndStdDF(df.groupby(COL_A))
    self.assertEqual(set(df_avg_std[cn.AVG].values.flatten()),
        set([0.5, 2.5]))

  def testMakeAvgAndStdDF2(self):
    if IGNORE_TEST:
      return
    DATA = [0.0236, 0.01511, 0.021991, 0.022334]
    df = pd.DataFrame({cn.RATE: DATA})
    key_name = 'HA2.152.08.03.D.CI'
    df[cn.KEY_ISOLATE] = key_name
    df_avg_std = ir.makeAvgAndStdDF(df.groupby(cn.KEY_ISOLATE))
    self.assertTrue(
        np.isclose(df_avg_std.loc[key_name, cn.AVG], 0.02075875))
    self.assertTrue(
        abs(df_avg_std.loc[key_name, cn.STD] - 0.001914) < 1e-5) 
    

class TestIsolateRegression(unittest.TestCase):

  def setUp(self):
    ir.IsolateRegression.resetDFS()
    self.regression = ir.IsolateRegression(
        ISOLATES, depvar=DEPVAR)
    self.cls = self.regression.__class__

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(self.regression._depvar, DEPVAR)
    self.assertTrue(isinstance(self.regression.dfs_coculture[cn.RATE],
        pd.DataFrame))
    self.assertTrue(isinstance(
        self.regression.dfs_coculture[cn.RATE], pd.DataFrame))
    self.assertEqual(
        [self.regression._isolate_dvh, self.regression._isolate_mmp], 
        ISOLATES)

  def testGetData(self):
    if IGNORE_TEST:
      return
    ir.IsolateRegression._getData()
    for depvar in [cn.RATE, cn.YIELD]:
      df = ir.IsolateRegression.dfs_coculture[depvar]
      self.assertTrue(helpers.isValidDataFrame(df,
          [cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP,
           cn.KEY_CULTURE, cn.DEPVAR]))

  def testKeyValues(self):
    if IGNORE_TEST:
      return
    isolate_pairs = self.regression.__class__.getKeyValues()
    trues = [Isolate.isSpecies(p[0], cn.SPECIES_MIX_DVH) and
             Isolate.isSpecies(p[1], cn.SPECIES_MIX_MMP)
             for p in isolate_pairs]
    self.assertTrue(all(trues))

  def testKeyValuesWithIsPermittedIsolate(self):
    if IGNORE_TEST:
      return
    isPermittedRow = lambda r: ISOLATES == [
        r[cn.KEY_ISOLATE_DVH], r[cn.KEY_ISOLATE_MMP]]
    isolate_pairs = self.regression.__class__.getKeyValues(
        isPermittedRow=isPermittedRow)
    self.assertEqual(len(isolate_pairs), 1)

  def testGetDepvars(self):
    if IGNORE_TEST:
      return
    df_full = self.regression._getDepvars(
        isIndex=lambda x: True)
    df_odd = self.regression._getDepvars(
        isIndex=lambda x: x % 2 == 1)
    df_even = self.regression._getDepvars(
        isIndex=lambda x: x % 2 == 0)
    union_cultures = set(df_odd[cn.KEY_CULTURE]).union(
        df_even[cn.KEY_CULTURE])
    full_cultures = set(df_full[cn.KEY_CULTURE])
    self.assertEqual(union_cultures, full_cultures)
    

class TestParametricIsolateRegression(unittest.TestCase):

  def setUp(self):
    ir.ParametricIsolateRegression.resetDFS()
    self.regression = ir.ParametricIsolateRegression(
        ISOLATES, depvar=DEPVAR)
    self.cls = self.regression.__class__

  def testEstimate(self):
    if IGNORE_TEST:
      return
    df = self.regression.estimate(isIndex=lambda x: True)
    self.assertTrue(helpers.isValidDataFrame(df,
        [rm.ALPHA_AVG, rm.ALPHA_STD, cn.KEY_ISOLATE_DVH,
         rm.BETA_AVG, rm.BETA_STD, cn.KEY_ISOLATE_MMP,
         rm.GAMMA_AVG, rm.GAMMA_STD],
        nan_columns=[rm.ALPHA_STD, rm.BETA_STD,
        rm.GAMMA_STD]))

  def testPredict(self):
    if IGNORE_TEST:
      return
    df = self.regression.estimate(isIndex=lambda x: x % 4 != 0)
    df = self.regression.predict(isIndex=lambda x: x % 4 == 0)
    self.assertTrue(helpers.isValidDataFrame(df,
        [cn.ESTIMATE, cn.OBSERVED, cn.KEY_CULTURE,
         cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP]))

  def testPlotAlphaEstimate(self):
    if IGNORE_TEST:
      return
    # smoke test only
    self.cls.plotAlphaEstimate(cn.YIELD, is_test=True)
    self.cls.plotAlphaEstimate(cn.RATE, is_test=True)
    

class TestNonParametricIsolateRegression(unittest.TestCase):

  def setUp(self):
    if IGNORE_TEST:
      return
    ir.NonParametricIsolateRegression.resetDFS()
    self.regression = ir.NonParametricIsolateRegression(
        ISOLATES, depvar=DEPVAR)
    self.cls = self.regression.__class__

  def testEstimate(self):
    if IGNORE_TEST:
      return
    df_estimate = self.regression.estimate(
        isIndex=lambda x: x % 3 != 0)
    self.assertTrue(helpers.isValidDataFrame(df_estimate,
        [cn.DEPVAR, cn.ESTIMATE, cn.KEY_ISOLATE_DVH, 
         cn.KEY_ISOLATE_MMP]))

  def testMakeEstimateDFS(self):
    if IGNORE_TEST:
      return
    cross = self.cls.makeEstimateDFS()
    self.assertEqual(len(cross), SIZE_OF_CROSS_RESULT)

  def testGetSmallResidualCultures(self):
    if IGNORE_TEST:
      return
    cultures10 = self.cls.getSmallResidualCultures(max_std=10)
    cultures3 = self.cls.getSmallResidualCultures(max_std=3)
    self.assertGreater(len(cultures10), len(cultures3))
    cultures2 = self.cls.getSmallResidualCultures(max_std=2)
    self.assertGreater(len(cultures3), len(cultures2))
    

class TestAncestralPairingIsolateRegression(unittest.TestCase):

  def setUp(self):
    ir.AncestralPairingIsolateRegression.resetDFS()
    self.regress = ir.AncestralPairingIsolateRegression(
        ISOLATES, depvar=DEPVAR, line=LINE)
    self.cls = self.regress.__class__

  def testGetFullIsPermittedRow(self):
    if IGNORE_TEST:
      return
    # Tests that isolate pairs are removed correctly
    isolate_pairs = util.removeDuplicatesFromList([
        (r[cn.KEY_ISOLATE_DVH], r[cn.KEY_ISOLATE_MMP])
        for _,r in self.cls.dfs_coculture[cn.RATE].iterrows()
        if r[cn.LINE] == LINE])
    #
    def test(isolate_pair, leave_out, size_adjustment):
      regress = ir.AncestralPairingIsolateRegression(
          isolate_pair, depvar=DEPVAR, line=LINE, leave_out=leave_out)
      constraint = regress._getFullIsPermittedRow()
      constrained_isolate_pairs = [
          (r[cn.KEY_ISOLATE_DVH], r[cn.KEY_ISOLATE_MMP])
          for _,r in self.cls.dfs_coculture[cn.RATE].iterrows() 
          if constraint(r)
          ]
      constrained_isolate_pairs = util.removeDuplicatesFromList(
          constrained_isolate_pairs)
      b = len(isolate_pairs) - size_adjustment == len(constrained_isolate_pairs)
      if not b:
        import pdb; pdb.set_trace()
      self.assertTrue(b)
    #
    size_dict = {cn.KEY_ISOLATE: 1, cn.KEY_CULTURE: 0}
    for isolate_pair in [ISOLATES, ISOLATES2]:
      for leave_out in [cn.KEY_ISOLATE, cn.KEY_CULTURE]:
        test(isolate_pair, leave_out, size_dict[leave_out])
    
  def testEstimate(self):
    if IGNORE_TEST:
      return
    def test(depvar, leave_out):
      self.regress = ir.AncestralPairingIsolateRegression(
          ISOLATES, depvar=depvar, line=LINE, 
          leave_out=leave_out)
      isIndex = lambda v: v != 0
      df = self.regress.estimate(isIndex=isIndex)
      self.assertTrue(helpers.isValidDataFrame(df,
          [cn.LINE, cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP,
           rm.ALPHA_AVG, rm.ALPHA_STD,
           rm.BETA_AVG, rm.BETA_STD, ir.COEF_DVH, ir.COEF_MMP,
           rm.GAMMA, cn.ESTIMATE, cn.RSQ], 
           nan_columns=[rm.ALPHA_STD, rm.BETA_STD]))
    #
    for depvar in [cn.RATE, cn.YIELD]:
      for leave_out in [cn.KEY_CULTURE, cn.KEY_ISOLATE]:
        test(depvar, leave_out)

  def testGetKeyValues(self):
    if IGNORE_TEST:
      return
    def test(line):
      key_values = self.cls.getKeyValues(line=line)
      self.assertTrue(len(key_values) > 0)
    #
    for line in ['HA2','HR2', 'UE3']:
      test(line)

  def testMakeEstimateDFS(self):
    if IGNORE_TEST:
      return
    cross = self.cls.makeEstimateDFS(line=cn.LINE_HA2)
    self.assertEqual(len(cross), SIZE_OF_CROSS_RESULT)
    df_ancestral = self.cls.dfs_ancestral[cn.RATE]
    #
    def test(line, leave_out):
      cross = self.cls.makeEstimateDFS(depvar=DEPVAR, line=line,
          leave_out=leave_out)
      if leave_out == cn.KEY_CULTURE:
        self.assertGreater(cross[cn.RSQ], 0.0)
      self.assertEqual(len(cross), SIZE_OF_CROSS_RESULT)
      self.assertTrue(helpers.isValidDataFrame(df_ancestral,
          [cn.AVG, cn.STD],
          nan_columns=[cn.STD]))
      self.assertFalse('index' in df_ancestral.columns)
      self.assertTrue(df_ancestral.index.name, cn.KEY_ISOLATE)
    #
    for line in ['HA2','HR2', 'UE3']:
      #for leave_out in [cn.KEY_CULTURE, cn.KEY_ISOLATE]:
      for leave_out in [cn.KEY_ISOLATE]:
        test(line, leave_out)

  def testIntegeration(self):
    if IGNORE_TEST:
      return
    cls_ir = self.cls
    depvar = cn.RATE
    max_std = 3.0
    isolate_pair = ("HA2.152.08.03.D.CI", "HA2.152.08.03.M.CI")
    leave_out_isolates = [isolate_pair]
    #
    def keepIsolates(row):
      if leave_out_isolates is None:
        return True
      if (row[cn.KEY_ISOLATE_DVH], row[cn.KEY_ISOLATE_MMP]) in leave_out_isolates:
        return False
      return True
    #   
    # Filter based on standard deviation of residuals
    cultures = cls_ir.getSmallResidualCultures(max_std)
    iprs = [lambda r: (r[cn.KEY_CULTURE] in cultures), 
        lambda r: keepIsolates(r) and (r[cn.KEY_CULTURE] in cultures),
        ]
    dfs_list = []
    for ipr in iprs:
      dfs_list.append(cls_ir.makeEstimateDFS(depvar=depvar, 
          isPermittedRow=ipr, 
          leave_out=cn.KEY_CULTURE, line='HA2'))
    self.assertEqual(len(dfs_list[0][cn.AVG]), 
        len(dfs_list[1][cn.AVG]) + 1)


if __name__ == '__main__':
    unittest.main()
