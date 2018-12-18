import __init__
import constants as cn
import helpers
import util
import binary_tree_classification as btc
import binary_tree_model as btm
from helpers_test_model import makeTree
import util_data as ud

import copy
import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
COL_A = 'a'
COL_B = 'b'
COL_C = 'c'
COL_D = 'd'
SIZE = 8
DF_X = pd.DataFrame({
    COL_A: [0, 1, 0, 1, 0, 1, 0, 1],
    COL_B: [0, 0, 1, 1, 0, 0, 1, 1],
    COL_C: [0, 0, 0, 0, 1, 1, 1, 1],
    })
DF_Y = pd.DataFrame({
    COL_A: [0, 1, 0, 1, 0, 1, 0, 1],
    })

################## TEST FUNCTIONS #################
class TestAll(unittest.TestCase):

  def setUp(self):
    self.cls = btc.BinaryTreeClassification
    self.bt_classify = self.cls(min_samples=1,
        is_validate=True)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(np.isnan(self.bt_classify.classify_value))

  def testCalcSplitAccuracy(self):
    if IGNORE_TEST:
      return
    ser = self.cls._calcSplitAccuracy(DF_X, DF_Y)
    self.assertEqual(set(ser.index), set([COL_A, COL_B, COL_C]))
    self.assertEqual(ser[COL_A], 1.0)
    for col in [COL_B, COL_C]:
      self.assertEqual(ser[COL_C], 0.5)

  def testFindSplitColumn(self):
    if IGNORE_TEST:
      return
    self.bt_classify.df_X = DF_X
    self.bt_classify.df_y = DF_Y
    col, score = self.bt_classify._findSplitColumn()
    self.assertEqual(col, COL_A)
    self.assertEqual(score, 1.0)
    #
    self.bt_classify.df_X = DF_X.applymap(lambda v: 1 - v)
    col, score = self.bt_classify._findSplitColumn()
    self.assertEqual(col, COL_A)
    self.assertEqual(score, 1.0)

  def testFitBasic(self):
    if IGNORE_TEST:
      return
    bt_classify = self.cls(min_samples=1, min_incr_score=0)
    bt_classify.fit(DF_X, DF_Y)
    for val in [0, 1]:
      self.assertEqual(bt_classify.trees[val].classify_value, 
          val)
      self.assertNotEqual(bt_classify.trees[val].split_column,
          cn.LEAF_NODE)
    self.assertEqual(len(bt_classify.findLeaves()),
        2**len(DF_X.columns))

  def testFitAndScoreScale(self):
    if IGNORE_TEST:
      return
    def getScore(predictor_probs, num_rows, noise_prob):
      probs = [0.2, 0.3, 0.5]
      probs = [1.0]
      bt_classify = self.cls(min_incr_score = 0.05)
      df_X, df_y, score = ud.makeNoisyOr(predictor_probs, num_rows, noise_prob)
      bt_classify.fit(df_X, df_y)
      classify_score = bt_classify.score(None, None)
      return classify_score
    #
    num_rows = 1000
    score1 = getScore([1.0], num_rows, 0.2)
    score2 = getScore([0.5, 0.5], num_rows, 0.2)
    score3 = getScore(np.repeat(0.2, 5), num_rows, 0.2)
    self.assertGreater(score1, score2)
    self.assertGreater(score2, score3)
    self.assertGreater(score3, 0)
    self.assertGreater(1, score1)

  def testScaleForNoisyBin(self):
    if IGNORE_TEST:
      return
    def getScore(num_rows, noise_prob):
      bt_classify = self.cls(max_depth=100,
          min_samples=1)
      df_X, df_y = ud.makeNoisyBin(num_rows, noise_prob)
      bt_classify.fit(df_X, df_y)
      classify_score = bt_classify.score(None, None)
      if False:
        print classify_score
        print bt_classify
        import pdb; pdb.set_trace()
      if np.isclose(noise_prob, 0):
        self.assertTrue(np.isclose(classify_score, 1.0))
      return classify_score
    #
    num_rows = 100
    SIZE = 5
    for _ in range(SIZE):
      score = getScore(num_rows, 0)
    for idx in range(SIZE):
      noise_prob = (1.0*idx) / SIZE
      new_score = getScore(num_rows, noise_prob)
      self.assertLess(new_score, score*1.01)
      score = new_score

  def testPredict(self):
    if IGNORE_TEST:
      return
    self.bt_classify.fit(DF_X, DF_Y)
    df = self.bt_classify.predict(DF_X)
    self.assertTrue(helpers.isValidDataFrame(df,
        [cn.AVG, cn.STD], nan_columns=[cn.STD]))
    trues = [v in [0.0, 1.0] for v in df[cn.AVG]]
    self.assertTrue(all(trues))

  def testGetParameterDF(self):
    if IGNORE_TEST:
      return
    self.bt_classify.fit(DF_X, DF_Y)
    df = self.bt_classify.getParameterDF()
    self.assertTrue(helpers.isValidDataFrame(df,
        DF_X.columns))

  def testPrune(self):
    if IGNORE_TEST:
      return
    """
    Test cases where should and should not prune
    Show get accuracy of 1.0 on NoisyBin if noise = 0.0.
    Create the tree
      a (root)
       - ab0 (node0)
         - leaf (leaf00)
         + leaf (leaf01)
       + leaf (leaf1)
    """
    ROOT = 'root'
    NODE0 = 'node0'
    LEAF00 = 'leaf00'
    LEAF01 = 'leaf01'
    LEAF01 = 'leaf01'
    LEAF1 = 'leaf1'
    def makeCustomTree(score_dict):
      root = makeTree(cls=btc.BinaryTreeClassification)
      node_dict = {ROOT: root}
      node_dict[LEAF1] = root.trees[1]
      node_dict[LEAF1].deleteSubtree()
      node0 = root.trees[0]
      node_dict[NODE0] = node0
      node_dict[LEAF00] = node0.trees[0]
      node_dict[LEAF00].deleteSubtree()
      node_dict[LEAF01] = node0.trees[1]
      node_dict[LEAF01].deleteSubtree()
      node_dict[NODE0].trees[1].deleteSubtree()
      for key in score_dict.keys():
        node_dict[key].node_score = score_dict[key]
      return node_dict
    def test(leaf00, leaf01, exp_score, exp_leaves, leaf1=0.8):
      score_dict = {
          ROOT: 0.6,
          NODE0: 0.7,
          LEAF00: leaf00,
          LEAF01: leaf01,
          LEAF1:  leaf1,
          }
      node_dict = makeCustomTree(score_dict)
      root = copy.deepcopy(node_dict[ROOT])
      root._prune()
      self.assertEqual(root.getScore(), exp_score)
      self.assertEqual(len(root.findLeaves()), exp_leaves)
    # Case 1: Nothing to prune
    test(1, 1, 0.9, 3)
    # Case 2: Prune LEAF00, LEAF01
    test(0.4, 0.5, 0.75, 2)
    # Case 3: All leaves
    test(0.4, 0.5, 0.6, 1, leaf1=0.5)

  def testHyperParameters(self):
    if IGNORE_TEST:
      return
    NUM_ROWS = 100
    NOISE_PROB = 0.0
    def getScore(min_samples=1, max_depth=10, min_incr_score=0):
      bt_classify = self.cls(
          min_samples=min_samples,
          max_depth=max_depth,
          min_incr_score=min_incr_score)
      df_X, df_y = ud.makeNoisyBin(NUM_ROWS, NOISE_PROB)
      bt_classify.fit(df_X, df_y)
      classify_score = bt_classify.getScore()
      return classify_score
    #
    for _ in range(5):
      self.assertLessEqual(getScore(min_incr_score=0.05),
          getScore(min_incr_score=0.0))
      self.assertLessEqual(getScore(min_samples=3), 
          getScore(min_samples=1))
      self.assertLessEqual(getScore(max_depth=3), 
          getScore(max_depth=10))
      

if __name__ == '__main__':
    unittest.main()
