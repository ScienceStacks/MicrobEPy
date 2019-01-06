"""Reorders rows and/or columns to group together categoricals."""

import util

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
import warnings


CLASS = "class"
SCORE = "score"

class DataframeSorter(object):

  def __init__(self, df):
    """
    :param pd.DataFrame df:
    """
    self._df = df

  def orderRows(self, max_clusters=10, is_convert_zero_to_nan=True):
    """
    Orders the rows in the dataframe so that:
      1. Like rows are together
      2. Rows with larger total values are towards the bottom
    :param int max_clusters: maximum number of clusters formed
    :param bool is_convert_zero_to_nan:
    :return pd.DataFrame: re-orders rows.
    """
    n_clusters = min(len(self._df), max_clusters)
    MULT = 0.001
    #
    df_result = self._df.copy()
    df_result = df_result.applymap(lambda v: 0 if np.isnan(v) else v)
    # Construct cluster groups
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      kmeans = KMeans(n_clusters=n_clusters).fit(df_result)
    df_result[CLASS] = kmeans.labels_
    # Compute group score
    df_result[SCORE] = df_result.sum(axis=1)
    dfg = pd.DataFrame(df_result.groupby(CLASS).sum())
    for idx, _ in dfg.iterrows():
      value = dfg.loc[idx, SCORE]
      sel = [s == idx for s in df_result[CLASS]]
      df_result.loc[sel, CLASS] = value
    df_result[CLASS] = df_result[CLASS] +   \
        MULT*df_result[SCORE]
    # Order the columns and reset zero values
    df_result = df_result.sort_values(CLASS)
    if is_convert_zero_to_nan:
      df_result = df_result.applymap(lambda v: np.nan if v==0 else v)
    # Cleanup the dataframe
    del df_result[CLASS]
    del df_result[SCORE]
    self._df = df_result
    #
    return self._df

  def orderColumns(self, **kwargs):
    """
    return pd.DataFrame
    """
    column_sorter = DataframeSorter(self._df.transpose())
    self._df = column_sorter.orderRows(**kwargs)
    self._df = self._df.transpose()
    return self._df

  def orderBoth(self, **kwargs):
    """
    Sorts by rows. Orders the columns in the same way.
    Assumes that the matrix has the same rows and columns.
    """
    self.orderRows(**kwargs)
    row_names = self._df.index.tolist()
    self._df = self._df[row_names]
    return self._df

  def deleteNanRowsAndColumns(self):
    """
    Deletes columns and rows that are all Nan.
    :return pd.DataFrame:
    """
    # Eliminate the NaN columns
    for col in self._df:
      length = len([v for v in self._df[col] if np.isnan(v)])
      if length == len(self._df):
        del self._df[col]
    # Eliminate the NaN rows
    for name in self._df.index:
      length = len([v for v in self._df.loc[name] if np.isnan(v)])
      if length == len(self._df.columns):
        self._df = self._df.drop([name])
    #
    return self._df
