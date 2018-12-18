import __init__
import helpers
from group_splitter import GroupSplitter
import constants as cn
import util

import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
COL_A = 'a'
COL_B = 'b'
COL_C = 'c'
COL_Y = 'y'
DATA = range(10)
DF_X = pd.DataFrame({
    COL_A: DATA,
    COL_B: DATA,
  })
DF_Y = pd.DataFrame({COL_Y: DATA})
DF_GROUP = pd.DataFrame({
    COL_A: ["a" + str(n) for n,_ in enumerate(DATA)],
    COL_B: ["b" + str(n) for n,_ in enumerate(DATA)],
  })
NUM_FOLDS = 2
  

class TestGroupSplitter(unittest.TestCase):

  def setUp(self):
    self.group_splitter = GroupSplitter(DF_X, DF_Y, DF_GROUP,
        num_folds=NUM_FOLDS)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    df = self.group_splitter._df_fold
    count1 = len(df[df[cn.FOLD] == 1])
    count0 = len(df[df[cn.FOLD] == 0])
    self.assertEqual(count0, count1)

  def testNext(self):
    if IGNORE_TEST:
      return
    train_y = []
    for dfs in self.group_splitter:
      self.assertEqual(
          set(dfs[cn.TRAIN_Y][COL_Y]).intersection(train_y), 
          cn.NULL_SET)
      train_y = dfs[cn.TRAIN_Y][COL_Y]
      for df in dfs.values():
        self.assertTrue(helpers.isValidDataFrame(
            df, df.columns))

  def testNextWithNone(self):
    if IGNORE_TEST:
      return
    self.group_splitter = GroupSplitter(DF_X, DF_Y, None)
    test_y = []
    for dfs in self.group_splitter:
      self.assertEqual(
          set(dfs[cn.TEST_Y][COL_Y]).intersection(test_y), 
          cn.NULL_SET)
      test_y.extend(dfs[cn.TEST_Y][COL_Y].tolist())
      for df in dfs.values():
        self.assertTrue(helpers.isValidDataFrame(
            df, df.columns))


if __name__ == '__main__':
    unittest.main()
