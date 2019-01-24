from microbepy.common import util
from microbepy.statistics.mutation_differential  \
    import MutationDifferential
from microbepy.common import constants as cn
from microbepy.common import helpers
from microbepy.plot.util_plot import PlotParms

import numpy as np
import os
import pandas as pd
import unittest

IGNORE_TEST = False

########################################
class TestMutationDifferential(unittest.TestCase):

  def setUp(self):
    self.differential = MutationDifferential(cn.RATE, cn.GGENE_ID)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(isinstance(self.differential.df_X, pd.DataFrame))
    self.assertGreater(len(self.differential.df_X), 10)
    #
    constraints = [lambda r: r[cn.LINE] == cn.LINE_HA2]
    differential = MutationDifferential(cn.RATE, cn.GGENE_ID,
        constraints=constraints)
    self.assertGreater(len(self.differential.df_X),
        len(differential.df_X))

  def testMakeCounts(self):
    if IGNORE_TEST:
      return
    def test(ser):
      self.assertTrue(isinstance(ser, pd.Series))
      trues = [(v >= 0) and (v <= len(self.differential.df_X))
          for v in ser]
      self.assertTrue(all(trues))
      tot = self.differential._tot_low + self.differential._tot_high
      self.assertEqual(tot, len(self.differential.df_X))
    #
    test(self.differential._ser_low)
    test(self.differential._ser_high)

  def testCalcSL(self):
    if IGNORE_TEST:
      return
    def testValid(sl):
      self.assertGreaterEqual(sl, 0)
      self.assertLessEqual(sl, 1.0)
    def testLess(pairs):
      sls = []
      for pair in pairs:
        sls.append(self.differential._calcSL(pair[0], pair[1]))
      self.assertLess(sls[0], sls[1])
    #
    sl1 = self.differential._calcSL(3, 0)
    sl2 = self.differential._calcSL(4, 0)
    self.assertGreater(sl1, sl2)
    #
    sl1 = self.differential._calcSL(15, 0)
    sl2 = self.differential._calcSL(15, 11)
    self.assertLess(sl1, sl2)
    #
    testValid(self.differential._calcSL(2, 2))
    testValid(self.differential._calcSL(0, 1))
    testValid(self.differential._calcSL(1, 0))
    #
    testLess([ [0, 10], [5, 5] ])
    testLess([ [2, 8], [5, 5] ])

  def testCalcProb(self):
    if IGNORE_TEST:
      return
    def test(tot_low, tot_high, num):
      probs = {}
      for nn in range(num+1):
        num_low = nn
        num_high = num - nn
        probs[nn] = MutationDifferential._calcProb(
            num_low, num_high, tot_low, tot_high)
      if not np.isclose(sum(probs.values()), 1.0):
        import pdb; pdb.set_trace()
      self.assertTrue(np.isclose(sum(probs.values()), 1.0))
    #
    for _ in range(10):
      tot_low = np.random.randint(5, 15)
      tot_high = np.random.randint(5, 15)
      num = np.random.randint(1, min(tot_low, tot_high))
      test(tot_low, tot_high, num)

  def testMakeDF(self):
    if IGNORE_TEST:
      return
    df1 = self.differential.makeDF()
    self.assertTrue(helpers.isValidDataFrame(df1,
        [cn.COUNT1, cn.COUNT2, cn.SIGLVL, cn.VALUE]))
    #
    differential = MutationDifferential(cn.RATE, cn.GGENE_ID,
        constraints=[lambda r: r[cn.LINE] == cn.LINE_UE3])
    df2 = differential.makeDF()
    self.assertGreater(len(df1), len(df2))

  def testMakeDFMedian(self):
    if IGNORE_TEST:
      return
    line = cn.LINE_HR2
    differential = MutationDifferential(cn.RATE, cn.GGENE_ID,
        constraints=[lambda r: r[cn.LINE] == line],
        is_median=False)
    df1 = differential.makeDF()
    differential = MutationDifferential(cn.RATE, cn.GGENE_ID,
        constraints=[lambda r: r[cn.LINE] == line],
        is_median=True)
    df2 = differential.makeDF()
    self.assertFalse(df1.equals(df2))

  def testPlot(self):
    # Smoke test
    line = cn.LINE_HA2
    depvar = cn.YIELD
    constraints = [lambda r: r[cn.LINE] == line]
    differential = MutationDifferential(depvar, cn.GGENE_ID,
        constraints=constraints,
        is_median=True,
        is_standardize_by_line = True,
        )
    parms = PlotParms()
    parms[cn.PLT_TITLE] = "%s, %s" % (line, depvar)
    differential.scatterKnob(parms=parms, is_plot=False, max_sl=0.95)

    
if __name__ == '__main__':
    unittest.main()
