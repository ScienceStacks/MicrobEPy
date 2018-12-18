import __init__
import constants as cn
import helpers
import util
import binary_tree_regression as btr
import util_data as ud

import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
COL_A = 'a'
COL_B = 'b'
COL_C = 'c'
SIZE_A = 3
SIZE_B = 2
DF_X = pd.DataFrame({
    COL_A: [1, 1, 1, 0],
    COL_B: [1, 1, 0, 0],
    })
DF_Y = pd.DataFrame({
    cn.VALUE: range(4)
    })
SMALL = 0.33
LARGE = 0.66
NAMES = ['a', 'b', 'c', 'd', 'e', 'f']
DEPTH = len(NAMES) - 1


################## TEST FUNCTIONS #################
class TestAll(unittest.TestCase):

  def setUp(self):
    self.cls = btr.BinaryTreeRegression
    self.bt_regress = btr.BinaryTreeRegression(min_samples=0,
        is_validate=True)

  def testCalcMSE(self):
    if IGNORE_TEST:
      return
    df = self.cls.calcMSE(DF_X, DF_Y)
    self.assertEqual(df.loc[COL_A],
        np.var(range(SIZE_A))*SIZE_A)
    self.assertEqual(df.loc[COL_B],
        np.var(range(SIZE_B))*SIZE_B)

  def testFindSplitColumn(self):
    if IGNORE_TEST:
      return
    column, mse = self.cls._findSplitColumn(DF_X, DF_Y)
    self.assertEqual(column, COL_B)
    self.assertEqual(mse, 1.0)

  def testFit1(self):
    if IGNORE_TEST:
      return
    self.bt_regress.fit(DF_X, DF_Y)
    for key in self.bt_regress.trees.keys():
      self.assertIsNotNone(self.bt_regress.trees[key])
    self.assertGreater(self.bt_regress.node_score, 
        self.bt_regress.split_score)

  def testFitScale(self):
    if IGNORE_TEST:
      return
    def makeStg(idx):
      return "X_%d" % idx
    #
    bt_regress = btr.BinaryTreeRegression(min_samples=10,
        is_validate=True)
    NUM_PREDICTORS = 5
    MAGNITUDE = 2
    means = [MAGNITUDE*v for v in range(1, NUM_PREDICTORS)]
    stds = np.repeat(1, NUM_PREDICTORS)
    num_rows = 1000
    num_extracols = 2
    def test():
      df_X, df_y = ud.makeGausianMixture(means, stds, num_rows, 
          num_extracols)
      self.bt_regress.fit(df_X, df_y)
      tree_str = str(self.bt_regress)
      # All of the predictor columns are present
      for mean in means:
        self.assertTrue(makeStg(mean) in tree_str)
      # Values appear in order of increasing value
      for idx in range(1, len(means)):
        pos1 = tree_str.index(makeStg(means[idx-1]))
        pos2 = tree_str.index(makeStg(means[idx]))
        self.assertGreater(pos1, pos2)
      self.assertGreater(self.bt_regress.score(df_X, df_y), 0.90)
    #
    for _ in range(5):
      test()

  def testScore(self):
    if IGNORE_TEST:
      return
    self.bt_regress.fit(DF_X, DF_Y)
    score = self.bt_regress.score(DF_X, DF_Y)
    self.assertTrue(isinstance(score, float))
      

if __name__ == '__main__':
    unittest.main()
