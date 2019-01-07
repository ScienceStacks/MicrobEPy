import microbepy_init
import constants as cn
import helpers
import util
from mutation_impact_plot import MutationImpactPlot

import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
IS_TEST = True
MUTATIONS = ['DVU1092', 'MMP1612']


################## TEST FUNCTIONS #################
class TestMutationImpactPlot(unittest.TestCase):

  def setUp(self):
    pass

  def testConsructor(self):
    if IGNORE_TEST:
      return
    mip = MutationImpactPlot(cn.RATE, cn.POSITION,
        is_test=IS_TEST)
    self.assertTrue(helpers.isValidDataFrame(mip.df_X,
        mip.df_X.columns))
    self.assertTrue(helpers.isValidDataFrame(mip.df_y,
        mip.df_y.columns))

  def testScatter(self):
    if IGNORE_TEST:
      return
    mip = MutationImpactPlot(cn.RATE, cn.GGENE_ID,
        is_test=IS_TEST)
    df_stat = mip.scatter(MUTATIONS)
    self.assertTrue(helpers.isValidDataFrame(df_stat,
        [cn.GROUP, cn.AVG, cn.STD, cn.COUNT]))
    #
    for col in [cn.GGENE_ID, cn.POSITION, cn.KEY_MUTATION]:
      mip = MutationImpactPlot(cn.RATE, col, is_test=IS_TEST)
    

if __name__ == '__main__':
    unittest.main()
