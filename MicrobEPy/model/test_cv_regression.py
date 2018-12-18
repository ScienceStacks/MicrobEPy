import __init__
import helpers
from cv_regression import CVLinearRegression, CVLassoRegression,  \
    CVRegression, CVForwardRegression, CVBinaryTreeRegression
import binary_tree_regression as btr
from group_splitter import GroupSplitter
import constants as cn
import util
import util_data as ud

import numpy as np
import pandas as pd
import random
from sklearn.model_selection import cross_val_predict, cross_val_score
from sklearn.metrics import r2_score
import unittest


IGNORE_TEST = False
COL_A = 'a'
COL_B = 'b'
COL_C = 'c'
COL_Y = 'y'
SIZE = 100
DATA_A = [random.normalvariate(0,1) for _ in range(SIZE)]
DATA_B = [random.normalvariate(0,1) for _ in range(SIZE)]
DF_X = pd.DataFrame({
    COL_A: DATA_A,
    COL_B: DATA_B,
  })
CONST_A = 7
CONST_B = 5
DF_Y = pd.DataFrame({
    COL_Y: [CONST_A*x + CONST_B*y for x,y in zip(DATA_A, DATA_B)]
    })
DF_GROUP = pd.DataFrame({
    COL_A: ["a" + str(n) for n,_ in enumerate(DATA_A)],
    COL_B: ["b" + str(n) for n,_ in enumerate(DATA_A)],
  })
NUM_FOLDS = 2


################### HELPERS ########################
def testCVR(cvr):
  """
  :param CVRegression cvr:
  :return dict: boolean values
  """
  result = {}
  cvr.fit()
  result['test3'] = helpers.isValidDataFrame(
      cvr.df_parameter, [cn.AVG, cn.STD, cn.COUNT])
  # df_parameter
  model = cvr.fitted_models[0]
  params = cvr.df_parameter.index.tolist()
  params.remove(cn.RSQ)
  for param in params:
    std = cvr.df_parameter.loc[param, cn.STD]
    result[param] = std < 0.01
  # df_predict
  result['test2'] = cvr.score > 0.95
  for key in result.keys():
    if not result[key]:
      import pdb; pdb.set_trace()
      pass
  return result
  

################### TEST CLASSES ########################
class TestCVRegression(unittest.TestCase):

  def setUp(self):
    self.g_splitter = GroupSplitter(DF_X, DF_Y, DF_GROUP,
        num_folds=NUM_FOLDS)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    cvr = CVLinearRegression(self.g_splitter)
    self.assertIsNone(cvr.df_parameter)

  def testFit(self):
    if IGNORE_TEST:
      return
    def test(cvr):
      results = testCVR(cvr)
      self.assertTrue(all(results.values()))
    #
    g_splitter = GroupSplitter(DF_X, DF_Y, DF_GROUP,
        num_folds=NUM_FOLDS)
    cvr = CVLassoRegression(g_splitter, alpha=0.001)
    test(cvr)
    #
    g_splitter = GroupSplitter(DF_X, DF_Y, DF_GROUP,
        num_folds=NUM_FOLDS)
    cvr = CVLinearRegression(g_splitter, fit_intercept=False, copy_X=True)
    test(cvr)


class TestCVForwardRegression(unittest.TestCase):

  def setUp(self):
    self.g_splitter = GroupSplitter(DF_X, DF_Y, DF_GROUP,
        num_folds=NUM_FOLDS)

  def testFit(self):
    if IGNORE_TEST:
      return
    df_X = DF_X.copy()
    df_X[COL_C] = [random.normalvariate(0,1) for _ in range(SIZE)]
    def test(cvr):
      results = testCVR(cvr)
      self.assertTrue(all(results.values()))
    #
    g_splitter = GroupSplitter(df_X, DF_Y, DF_GROUP,
        num_folds=NUM_FOLDS)
    cvr = CVForwardRegression(g_splitter,
        fit_intercept=False, copy_X=True)
    test(cvr)
    self.assertEqual(set(cvr.df_parameter.index),
        set([COL_A, COL_B, cn.RSQ]))


class TestCVBinaryTreeRegression(unittest.TestCase):

  def setUp(self):
    self.df_X, self.df_y = ud.makeGausianMixture([5, 10, 20],
        [1, 1, 1], 100, 0)
    self.df_group = ud.makeGroups(len(self.df_X))
    self.g_splitter = GroupSplitter(self.df_X, self.df_y, 
        self.df_group, num_folds=NUM_FOLDS)
    self.bt_regress = CVBinaryTreeRegression(self.g_splitter)

  def testPredictDF(self):
    if IGNORE_TEST:
      return
    _, self.bt_regress._fitted_model =  \
        self.bt_regress.fitDF(self.df_X, self.df_y)
    self.bt_regress.predictDF(self.df_X, self.df_y)
    df = self.bt_regress.predictDF(self.df_X, self.df_y)
    self.assertTrue(helpers.isValidDataFrame(df,
        [cn.OBSERVED, cn.ESTIMATE, cn.RESIDUAL]))

  def testFitDF(self):
    if IGNORE_TEST:
      return
    df, _ = self.bt_regress.fitDF(self.df_X, self.df_y)
    columns = self.df_X.columns.tolist()
    columns.append(cn.RSQ)
    self.assertTrue(helpers.isValidDataFrame(df, columns))

  def testFit(self):
    #
    def test(cvr):
      results = testCVR(cvr)
      self.assertTrue(all(results.values()))
    #
    cvr = CVBinaryTreeRegression(self.g_splitter)
    test(cvr)


if __name__ == '__main__':
    unittest.main()
