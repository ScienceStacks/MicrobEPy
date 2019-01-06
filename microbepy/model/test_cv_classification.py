import __init__
import helpers
from cv_classification import CVBinaryTreeClassification
import binary_tree_classification as btc
from group_splitter import GroupSplitter
import constants as cn
import util
import util_data as ud

import numpy as np
import pandas as pd
import random
import unittest


IGNORE_TEST = False
COL_A = 'a'
COL_B = 'b'
COL_C = 'c'
DF_X = pd.DataFrame({
    COL_A: [1, 1, 1, 0],
    COL_B: [1, 1, 0, 0],
    })
DF_Y = pd.DataFrame({
    cn.VALUE: [1, 1, 0, 0],
    })
DF_GROUP = pd.DataFrame({
    COL_A: ["a" + str(n) for n in range(len(DF_X))],
    COL_B: ["b" + str(n) for n in range(len(DF_X))],
  })
NUM_FOLDS = len(DF_X)


######################## HELPERS ####################
def calcScoreNoisyBin(num_rows, noise_prob, **kwargs):
  df_X, df_y = ud.makeNoisyBin(num_rows, noise_prob)
  g_splitter = GroupSplitter(df_X, df_y, None)
  if not 'min_incr_score' in kwargs.keys():
    kwargs['min_incr_score'] = 0.02
  cvc = CVBinaryTreeClassification(g_splitter, 
      **kwargs)
  cvc.fit()
  return cvc
  

################### TEST CLASSES ########################
class TestCVClassification(unittest.TestCase):

  def setUp(self):
    self.g_splitter = GroupSplitter(DF_X, DF_Y, DF_GROUP,
        num_folds=NUM_FOLDS)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    cvc = CVBinaryTreeClassification(self.g_splitter)
    self.assertIsNone(cvc.df_parameter)

  def testFit(self):
    if IGNORE_TEST:
      return
    #
    cvc = CVBinaryTreeClassification(self.g_splitter,
        max_node_score=1.0, min_samples=1)
    cvc.fit()
    self.assertEqual(cvc.df_parameter.loc[COL_B, cn.AVG], 1.0)
    self.assertTrue(helpers.isValidDataFrame(cvc.df_parameter,
        [cn.AVG, cn.COUNT, cn.STD]))
    self.assertTrue(helpers.isValidDataFrame(cvc.df_predict,
        [cn.ESTIMATE, cn.OBSERVED, cn.RESIDUAL]))
    self.assertEqual(len(cvc.fitted_models), NUM_FOLDS)

  def testPredictDF(self):
    if IGNORE_TEST:
      return
    cvc = CVBinaryTreeClassification(self.g_splitter,
        max_node_score=1.0, min_samples=1)
    cvc.fit()
    self.assertTrue(helpers.isValidDataFrame(cvc.df_predict,
        [cn.ESTIMATE, cn.OBSERVED, cn.RESIDUAL]))
  
  def testCalcPredictScore(self):
    if IGNORE_TEST:
      return
    cvc = CVBinaryTreeClassification(self.g_splitter,
        max_node_score=1.0, min_samples=1)
    cvc.fit()
    score = cvc.calcPredictScore()
    self.assertEqual(score, 1.0)

  def testScale(self):
    if IGNORE_TEST:
      return
    #
    def calcScore(predictor_probs, num_rows, noise_prob):
      df_X, df_y, score = ud.makeNoisyOr(predictor_probs, 
          num_rows, noise_prob)
      g_splitter = GroupSplitter(df_X, df_y, None)
      cvc = CVBinaryTreeClassification(g_splitter)
      cvc.fit()
      return cvc
    # Do better with either low or large predictor probability of
    # a dependent event
    num_rows = 100
    cvc = calcScore([1], num_rows, 0.0)
    self.assertTrue(np.isclose(cvc.score, 1))
    cvc1 = calcScore([1.0], num_rows, 0.2)
    cvc2 = calcScore([0.5, 0.5], num_rows, 0.2)
    self.assertGreater(cvc1.score, cvc2.score)

  def testScaleNoisyBin(self):
    if IGNORE_TEST:
      return
    num_rows = 50
    cvcs = []
    for prob in [0.5, 0.8]:
      cvcs.append(calcScoreNoisyBin(num_rows, prob))
    trues = [cvc.score < 0.8 for cvc in cvcs]
    self.assertTrue(all(trues))

  def testScaleNoisyBin2(self):
    if IGNORE_TEST:
      return
    num_rows = 50
    cvc = calcScoreNoisyBin(num_rows, 0.0,
        max_depth=10, min_samples=1, min_incr_score=0.0)
    self.assertEqual(cvc.df_parameter.loc[cn.ACC, cn.AVG], 1.0)
    self.assertEqual(len(cvc.fitted_models[0].findLeaves()), 
        num_rows-1)


if __name__ == '__main__':
    unittest.main()
