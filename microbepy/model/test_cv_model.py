from microbepy.common import constants as cn
from microbepy.common import helpers
from microbepy.common import util
from microbepy.data import util_data as ud
from microbepy.model.cv_regression import CVLinearRegression
from microbepy.model.group_splitter import GroupSplitter

import numpy as np
import pandas as pd
import random
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
def dotestCVR(cvr):
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
class TestCVModel(unittest.TestCase):

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
      results = dotestCVR(cvr)
      self.assertTrue(all(results.values()))
    #
    g_splitter = GroupSplitter(DF_X, DF_Y, DF_GROUP,
        num_folds=NUM_FOLDS)
    cvr = CVLinearRegression(g_splitter, fit_intercept=False, copy_X=True)
    test(cvr)

  def testFindSignificantParameters(self):
    if IGNORE_TEST:
      return
    g_splitter = GroupSplitter(DF_X, DF_Y, DF_GROUP,
        num_folds=NUM_FOLDS)
    cvr = CVLinearRegression(g_splitter, fit_intercept=False, copy_X=True)
    cvr.fit()
    parameters = cvr.findSignificantParameters()
    self.assertEqual(set(parameters), set([COL_A, COL_B]))
    


if __name__ == '__main__':
    unittest.main()
