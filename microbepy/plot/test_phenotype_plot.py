import __init__
import constants as cn
import helpers
import util
from phenotype_plot import PhenotypePlot

import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
IS_PLOT = False


################## TEST FUNCTIONS #################
class TestPhenotypePlot(unittest.TestCase):

  def setUp(self):
    self.plot = PhenotypePlot(is_plot=IS_PLOT, constraints=None)

  def testConsructor(self):
    if IGNORE_TEST:
      return
    self.assertGreater(len(self.plot.lines), 0)

  def testScatter(self):
    self.plot.scatter(is_errorbars=True)

  def testGetPoints(self):
    df = self.plot.getPoints()
    self.assertTrue(helpers.isValidDataFrame(df,
        [cn.RATE, cn.YIELD]))
    

if __name__ == '__main__':
    unittest.main()
