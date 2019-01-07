"""Splits Data into training and test data by groups."""

import microbepy_init
import constants as cn

import numpy as np
import pandas as pd

class GroupSplitter:
  """
  This class is an iterator. Usage:
    gsplitter = GroupSplitter(X, y)
    for dfs in gsplitter:
      <use dfs>
  """

  def __init__(self, df_X, df_y, df_group, num_folds=10):
    """
    :param pd.DataFrame df_X: predictor variables indexed by observation
    :param pd.DataFrame df_y: dependent variable indexed by observation
    :param pd.DataFrame df_group: columns define groups for fold,
                                  indexed by observation
                                  if none, use the index.
    :param int num_folds: number of folds (groups)
    """
    self._df_X = df_X
    self._df_y = df_y
    if df_group is None:
      self._df_group, num_folds = self._makeDefaultGroupDF()
    else:
      self._df_group = df_group
    self._num_folds = num_folds
    self._df_fold = self._makeFoldDF()
    self._fold = 0

  def _makeDefaultGroupDF(self):
    """
    Constructs a group DF with a key for each row.
    """
    num_folds = len(self._df_X)
    names = ["a%d" % n for n in range(num_folds)]
    df = pd.DataFrame({'a': names})
    return df, num_folds

  def _makeFoldDF(self):
    """
    Creates a dataframe that relates observation index to the fold.
    :return pd.DataFrame: cn.FOLD, indexed as self._df_X, _y
    """
    # Extract the columns that define rows in the same fold (a group)
    # Retain the mapping between the index and the group
    # Get lists of the fold groups
    columns = self._df_group.columns.tolist()
    groups = self._df_group.groupby(columns).groups
    num_groups = len(groups.keys())
    pairs = zip(groups.keys(), np.random.permutation(num_groups))
    index_fold_dict = {k: n % self._num_folds for k,n in pairs}
    # Associate each index with a fold
    group_ids = []
    for _, row in self._df_group.iterrows():
      if len(columns) == 1:
        key = row[columns[0]]
      else:
        key = tuple([row[c] for c in columns])
      group_ids.append(index_fold_dict[key])
    df_fold = pd.DataFrame({
      cn.FOLD: group_ids,
      })
    df_fold.index = self._df_group.index
    return df_fold

  def __iter__(self):
    """
    Returns the next split.
    :return dict-pd.DataFarme: TRAIN_X, TRAIN_Y, TEST_X, TEST_Y;
      all have their original indices
    """
    while self._fold < self._num_folds:
      dfs = {}
      sel_test = self._df_fold[
          self._df_fold[cn.FOLD] == self._fold].index.tolist()
      sel_train = self._df_fold[
          self._df_fold[cn.FOLD] != self._fold].index.tolist()
      dfs[cn.TEST_X] =  self._df_X.loc[sel_test].copy()
      dfs[cn.TEST_Y] =  self._df_y.loc[sel_test].copy()
      dfs[cn.TRAIN_X] =  self._df_X.loc[sel_train].copy()
      dfs[cn.TRAIN_Y] =  self._df_y.loc[sel_train].copy()
      self._fold += 1
      yield dfs
