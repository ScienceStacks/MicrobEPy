"""Computes Simple Statistics From the Data"""

import __init__
from api import Api
import constants as cn
import util

import pandas as pd
import numpy as np

class CoStatistics(object):

  def __init__(self, api_object=None):
    if api_object is None:
      api_object = Api()
    self.api = api_object

  def makeCountDF(self, 
      row_label=cn.KEY_MUTATION, 
      for_label=cn.KEY_ISOLATE,
      col_label=None,
      **kwargs):
    """
    Creates a DataFrame with the rows identifyed by the column row_label
    counting occurrences of the column for_label partitioning counts
    by distanct values of col_label.
    :param str row_label: name of the column to use as values for row labels
    :param str for_label: name of the column for which values are counted
    :param str col_label: name of the column whose values are used to
        partition the counts
    :return pd.DataFrame: columns: row_label, values in col_label
    """
    if col_label is None:
      df = self.api.makeDF(columns=[row_label, for_label], **kwargs)
      df_result = df.groupby([row_label]).count()
      df_result.rename(columns={for_label: cn.COUNT},
          inplace=True)
    else:
      df = self.api.makeDF(columns=[row_label, for_label, col_label], **kwargs)
      df_g = df.groupby([row_label, col_label]).count()
      df_g.reset_index(inplace=True)
      df_result = pd.DataFrame(df_g.pivot_table(index=row_label, values=for_label,
          columns=col_label, fill_value=0, aggfunc=np.sum))
    df_result.reset_index(inplace=True)
    util.changeColumnValues(df_result,
        lambda c, x: 0 if util.isNull(x) else x)
    util.cleanDF(df_result)
    return df_result
