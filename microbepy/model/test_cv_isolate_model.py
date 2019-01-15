from microbepy.common import helpers
from microbepy.common import util
from microbepy.common import constants as cn
from microbepy.model.cv_isolate_model import CVIsolateModel
from microbepy.model import isolate_model as rm
from microbepy.model import isolate_regression as ir

import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False

################ SUPPORTING CLASS #############
class IsolateModelTester(rm.IsolateModel):
  """
  Artificial class used for testing
  """
  COL_A = cn.KEY_ISOLATE_DVH
  COL_B = cn.KEY_ISOLATE_MMP
  COL_C = 'c'
  key_names = [COL_A, COL_B]
  key_values = {k: ['x', 'x', 'y', 'y'] for k in key_names}
  data = dict(key_values)
  data[COL_C] = range(len(key_values[key_names[0]]))
  df_data = pd.DataFrame(data)

  @classmethod
  def getKeyValues(cls, isPermittedRow=rm.ISPERMITTEDROW):
    """
    :return list-of-object: list of values of keys
    """
    result = []
    for idx in range(len(cls.key_values[cls.key_names[0]])):
      result.append(
          (cls.key_values[cls.key_names[0]][idx],
           cls.key_values[cls.key_names[1]][idx]))
    return result

  def getDataSize(self):
    """
    :return int: Number of data items (typically replications)
    """
    cls = self.__class__
    return len(cls.df_data.index)

  def estimate(self, isIndex=lambda x: True, **kwargs):
    """
    :param function isIndex:
       arg: int; returns bool
    :return pd.DataFrame
      Rows of original dataframe that are selected by isIndex
    """
    cls = self.__class__
    rows = []
    for idx, row in cls.df_data.iterrows():
      if isIndex(idx):
        rows.append(row)
    return pd.DataFrame(rows)

  def predict(self, isIndex=lambda x: True):
    """
    :param Function isIndex:
    :return pd.DataFrame:
      key columns
      cn.ESTIMATE - estimated value
      cn.OBSERVED
      Rows of original dataframe that are selected by isIndex
    """
    cls = self.__class__
    dfs = []
    for key_values in cls.df_data.groupby(
        cls.key_names).indices.keys():
      df = pd.DataFrame()
      df[cls.COL_A] = [key_values[0]]
      df[cls.COL_B] = [key_values[1]]
      dfs.append(df)
    result = pd.concat(dfs, sort=True)
    result[cn.ESTIMATE] = range(len(result.index))
    result[cn.OBSERVED] = result[cn.ESTIMATE]
    return result


################ CLASS UNDER TEST ############
class TestCVIsolateModel(unittest.TestCase):

  def setUp(self):
    self.regr_cls = IsolateModelTester
    self.cross = CVIsolateModel(self.regr_cls, 2)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(
        len(self.regr_cls.key_values[self.regr_cls.COL_A]),
        len(self.cross._key_values))

  def testEstimate(self):
    if IGNORE_TEST:
      return
    estimate = self.cross.estimate()
    self.assertEqual(set(estimate.keys()), 
        set([cn.AVG, cn.STD, cn.RSQ, cn.RESIDUAL]))
    self.assertEqual(estimate[cn.RSQ], 1.0)
    for key in [cn.AVG, cn.STD, cn.RESIDUAL]:
      self.assertTrue(helpers.isValidDataFrame(estimate[key],
          estimate[key].columns,
          nan_columns=estimate[key].columns))

  def testIntegration(self):
    if IGNORE_TEST:
      return
    cross = CVIsolateModel(ir.AncestralPairingIsolateRegression,
        leave_out=cn.KEY_ISOLATE, line='HA2')
    estimate = cross.estimate()
    for key in [cn.AVG, cn.RESIDUAL]:
      nan_columns = [rm.ALPHA_STD, rm.BETA_STD]
      self.assertTrue(helpers.isValidDataFrame(estimate[key],
          estimate[key].columns, nan_columns=nan_columns))
    rsqs = estimate[cn.AVG][cn.RSQ].tolist()
    self.assertEqual(len(rsqs), len(set(rsqs)))

  def testIntegration2(self):
    if IGNORE_TEST:
      return
    rsqs = {}
    for depvar in [cn.RATE, cn.YIELD]:
      dfs  = ir.NonParametricIsolateRegression.makeEstimateDFS(
          depvar=depvar)
      rsqs[depvar] = dfs[cn.RSQ]
    self.assertGreater(rsqs[cn.RATE], rsqs[cn.YIELD])

  def testIntegeration3(self):
    import matplotlib.pyplot as plt
    depvar = cn.YIELD
    line = cn.LINE_HA2
    for line in cn.LINE_CIS:
      cls = ir.AncestralPairingIsolateRegression
      cross = CVIsolateModel(cls, depvar=depvar,
          line=line, leave_out=cn.KEY_ISOLATE)
      result = cross.estimate()
      df = result[cn.RESIDUAL]
      plt.scatter(df[cn.OBSERVED], df[cn.ESTIMATE])
      plt.title("RSQ=%2.4f" % result[cn.RSQ])
      plt.xlabel('Observed')
      plt.ylabel('Estimated')
      #plt.show()
      

if __name__ == '__main__':
    unittest.main()
