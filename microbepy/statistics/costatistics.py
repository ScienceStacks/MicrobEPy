"""Computes Simple Statistics From the Data"""

import microbepy_init
import constants as cn
import util
import util_data

import pandas as pd
import numpy as np


class CoStatistics(object):

  def __init__(self, df_data=None):
    """
    :param pd.DataFrame df_data: genotype_phenotype data
    """
    if df_data is None:
      df_data = util_data.makeIsolateData(is_separate_species=False)
    self.df_data = df_data

  def makeCountDF(self, 
      row_label=cn.KEY_MUTATION, 
      for_label=cn.KEY_ISOLATE,
      col_label=None,
      constraints=None):
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
    df = util.selectRows(self.df_data, constraints)
    if col_label is None:
      df = df[[row_label, for_label]].copy()
      df_result = df.groupby([row_label]).count()
      df_result.rename(columns={for_label: cn.COUNT},
          inplace=True)
    else:
      df = df[[row_label, for_label, col_label]].copy()
      df_g = df.groupby([row_label, col_label]).count()
      df_g.reset_index(inplace=True)
      df_result = pd.DataFrame(df_g.pivot_table(index=row_label, values=for_label,
          columns=col_label, fill_value=0, aggfunc=np.sum))
    df_result.reset_index(inplace=True)
    util.changeColumnValues(df_result,
        lambda c, x: 0 if util.isNull(x) else x)
    util.cleanDF(df_result)
    return df_result
