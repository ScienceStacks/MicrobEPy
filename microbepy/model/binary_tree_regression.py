"""Regression tree for binary valued variables."""

from microbepy.common import constants as cn
from microbepy.common import util
from microbepy.model import binary_tree_model

import copy
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score

# Hyperparameters
MIN_RSQ = 0.05  # Minimum value of RSQ to extend a tree
#
NO_SPLIT_MIN_RSQ = "min_rsq"


class BinaryTreeRegression(binary_tree_model.BinaryTreeModel):

  def __init__(self, min_rsq=MIN_RSQ, **kwargs):
    """
    :param float min_rsq: minimal incremental improvement in RSQ
    :param dict kwargs: optional arguments passed to parent
    """
    super(BinaryTreeRegression, self).__init__(**kwargs)
    self._min_rsq = min_rsq

  def _printStats(self):
    col_y = self.df_y.columns[0]
    mse = np.var(self.df_y[col_y]) * len(self.df_y)
    stg = "MSE=%f AVG=%f STD=%f RSQ=%1.2f" % (
        mse, self.mean, self.std, self.incremental_score)
    return stg

  @staticmethod
  def calcMSE(df_X, df_y):
    """
    Calculates the MSE for the rows selected by each column of df_X.
    :param pd.DataFrame df_X:
    :param pd.DataFrame df_y:
    :return pd.DataFrame: cn.VALUE, row index is column
    """
    df_count = df_X.sum()
    width = len(df_X.columns)
    df_y_wide = pd.concat([df_y]*width, axis=1, sort=True)
    df_y_wide.columns = df_X.columns
    df_val = df_X * df_y_wide
    df_sum = df_val.sum()
    df_ssq = df_val * df_val
    df_ssq = df_ssq.sum()
    df_result = df_ssq - df_sum*df_sum/df_count
    return df_result

  @classmethod
  def _findSplitColumn(cls, df_X, df_y):
    """
    :param pd.DataFrame df_X: columns with binary values
    :param pd.DataFrame df_y: single column with float
    :return str, float: split column, mse
    """
    # Find the column that minimizes MSE
    df_X_c = df_X.applymap(lambda v: 1 - v)
    df_mse = cls.calcMSE(df_X, df_y) + cls.calcMSE(df_X_c, df_y)
    df_mse_sort = df_mse.sort_values()
    mse = df_mse_sort[0]
    split_column = df_mse_sort.index[0]
    #
    return split_column, mse

  def fit(self, df_X, df_y, is_first_call=True):
    """
    Recursively find splits to form tree.
    Updates instance variables as appropriate.
    :param pd.DataFrame df_X: columns with binary values
    :param pd.DataFrame df_y: single column with float
    :param bool is_first_call: indicates first call to fit
    Updates self.split_column, self.split_score, self.trees
    """
    self.validate()
    # Update state
    calcNodeScore = lambda: self.count*np.var(
        self.df_y[self.df_y.columns[0]])
    self._updateState(df_X, df_y, calcNodeScore)
    #
    self.split_column = cn.LEAF_NODE  # Assume leaf node
    # Find the splitable columns
    columns = self.getSplitColumns()
    if len(columns) == 0:
      self.no_split = binary_tree_model.NO_SPLIT_NO_COLUMNS
      return
    self.df_X = self.df_X[columns].copy()
    #
    split_column, split_score = self._findSplitColumn(
        self.df_X, self.df_y)
    if self.node_score > 0:
      rsq = 1 - split_score/self.node_score
    else:
      rsq = 0
    if rsq < self._min_rsq:
      self.no_split = NO_SPLIT_MIN_RSQ
      return
    # Acceptable split
    self.split_column = split_column
    self.split_score = split_score
    self.incremental_score = rsq
    # Recursively do splits
    self._doSplits()
    self.validate()

  def score(self, df_X, df_y):
    """
    Computes the R2 score
    :param pd.DataFrame df_X:
    :param pd.DataFrame df_y:
    :return float:
    """
    df_predict = self.predictDF(df_X, df_y)
    score = 1 - np.var(df_predict[cn.RESIDUAL])  \
        / np.var(df_predict[cn.OBSERVED])
    return score

  def calcPredictValue(self):
    """
    Provides the predicted value for the node.
    :return float:
    """
    return self.mean
