import microbepy_init
import helpers
import isolate_model as im
import util

import copy
import collections
import constants as cn
import numpy as np
import pandas as pd
import os
import unittest


IGNORE_TEST = False
DEPVAR = cn.RATE
    

class TestIsolateModel(unittest.TestCase):

  def testMakeAncestralDFS(self):
    if IGNORE_TEST:
      return
    dfs = im.IsolateModel._makeAncestralDFS()
    # Examine the dataframe
    for key in list(dfs.keys()):
      df = dfs[key]
      self.assertTrue(helpers.isValidDataFrame(df, [cn.AVG, cn.STD],
          nan_columns=[cn.STD]))
      self.assertTrue(df.index.name, cn.KEY_ISOLATE)
    #
    for key in list(dfs.keys()):
      df = dfs[key]
      self.assertTrue(helpers.isValidDataFrame(df,
          [cn.AVG, cn.STD], nan_columns=[cn.STD]))
    #
    key_name = 'HA2.152.08.03.D.CI'
    df = dfs[cn.RATE]
    self.assertTrue(
        abs(df.loc[key_name, cn.AVG] - 0.020774) < 1e-5) 
    self.assertTrue(
        abs(df.loc[key_name, cn.STD] - 0.001912) < 1e-5) 

  def testMakeCocultureDFS(self):
    if IGNORE_TEST:
      return
    dfs = im.IsolateModel._makeCocultureDFS()
    for key in list(dfs.keys()):
      df = dfs[key]
      self.assertTrue(helpers.isValidDataFrame(df,
          [cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP,
           cn.KEY_CULTURE, cn.DEPVAR]))

  def testMakeCultureIsolateDF(self):
    if IGNORE_TEST:
      return
    df = im.IsolateModel.makeCultureIsolateDF()
    self.assertTrue(helpers.isValidDataFrame(df,
        [cn.KEY_CULTURE, cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP]))



if __name__ == '__main__':
    unittest.main()
