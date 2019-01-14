from microbepy.common import constants as cn
from microbepy.common import util
from microbepy.common import helpers
from microbepy.common.mutation_context import MutationContext
from microbepy.search.mutation_combination import MutationCombination

import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
DEPVAR = cn.RATE
THIS_MUTATION_COLUMN = cn.POSITION
MUTATION_CONTEXT = MutationContext(DEPVAR, THIS_MUTATION_COLUMN)


class TestMutationCombination(unittest.TestCase):

  def setUp(self):
    self.combination = MutationCombination(MUTATION_CONTEXT)

  def testConstructor(self): 
    if IGNORE_TEST:
      return
    self.assertEqual(len(self.combination._constraints), 0)

  def testDoBasic(self):
    if IGNORE_TEST:
      return
    def test(size, is_resample=False, lines=None):
      df_min, df_max =   \
          self.combination.do(size, is_tstat=False, is_resample=is_resample,
          lines=lines)
      combinations = df_min[cn.MUTATIONS]
      trues = [len(c) in range(1,size+1) for c in combinations]
      self.assertTrue(all(trues))
      for df in [df_min, df_max]:
        self.assertTrue(helpers.isValidDataFrame(df,
            [cn.MUTATIONS, cn.SL_TSTAT,
            cn.SL_RESAMPLE, cn.LINE],
            nan_columns= [cn.SL_TSTAT, cn.SL_TSTAT,
            cn.SL_RESAMPLE]))
        self.assertFalse(all([np.isclose(x,y)
            for x, y in zip(df_min[cn.AVG], df_max[cn.AVG])]))
    #
    test(1, lines=["HA2", "HR2"])
    test(2, lines=["HA2"])


if __name__ == '__main__':
    unittest.main()
