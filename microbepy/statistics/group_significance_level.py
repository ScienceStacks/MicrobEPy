"""
Computes Significance Levels for the occurrence
of observations grouped by mutations being present or absent.
"""

from microbepy.common import constants as cn
from microbepy.common import util
from microbepy.data import util_data as ud

import numpy as np
import pandas as pd
from scipy.stats import t as t_dist


CDF_PROB_MIN = "cdf_prob_min"
CDF_PROB_MIN_ONE_MINUS = "cdf_prob_min_one_minus"
CDF_PROB_MAX = "cdf_prob_max"


######################################
# CLASSES
######################################
class SignificanceLevel(object):
  """  
  Information associated with a group significance level
  """  

  def __init__(self, sl_tstat=np.nan, sl_resample=np.nan,
      mutations=None, group=None, avg=np.nan, count=np.nan):
    self.avg = avg  # Average value of observations for the group
    self.group = group  # Binary representation of mutations present/absent
    self.count = count  # Count of items in the group
    self.mutations = mutations  # Ordered list of mutations for group
    self.sl_tstat = sl_tstat  # Significance level found by t-statistic
    self.sl_resample = sl_resample  # Significance level found by resample

  def __repr__(self):
    def assignValue(num):
      if np.isnan(num):
        return str(np.nan)
      else:
        return "%2.4f" % num
    sl_tstat = assignValue(self.sl_tstat)
    sl_resample = assignValue(self.sl_resample)
    avg = assignValue(self.avg)
    count = assignValue(self.count)
    stg = '''
    sl_tstat: %s, sl_resample: %s, group: %s\n
    mutations: %s\n
    avg: %s, count: %s
    ''' % (sl_tstat, sl_resample, self.group, str(self.mutations), avg, count)
    return stg

######################################
class GroupSignificanceLevel(object):

  def __init__(self, provider, mutations):
    """
    :param ModelDataProvider provider: data source
    :param list-str mutations: subset of mutations know to provider
    """
    self._provider = provider
    self._mutations = mutations
    self.df_group = ud.makeMutationGroupDF(
        self._provider.df_X, self._provider.df_y, 
        self._mutations)
    self.df_stat = self.makeGroupStatisticDF()
    self.sl_min = SignificanceLevel()
    self.sl_max = SignificanceLevel()
    self._initSignificanceLevels()

  def makeGroupStatisticDF(self,
      group_column=cn.GROUP, 
      value_column=cn.DEPVAR):
    """
    :param str group_column: column with the grouping data
    :param str value_column: column with the data values
    :return pd.DataFrame: cn.GROUP, cn.AVG, cn.STD, cn.COUNT
        The STD is the standard deviation of the mean
    """
    df_result = pd.DataFrame()
    df_result[cn.AVG] = self.df_group.groupby(
        group_column).mean()[value_column]
    df_result[cn.STD] = self.df_group.groupby(
        group_column).std()[value_column]
    df_result[cn.COUNT] = self.df_group.groupby(
        group_column).count()[value_column]
    df_result[cn.STD] = df_result[cn.STD]  \
        / df_result[cn.COUNT].apply(lambda v: np.sqrt(v))
    for col in [cn.STD]:
      df_result[col] = df_result[col].apply(
          lambda v: 0 if np.isnan(v) else v)
    df_result = df_result.reset_index()
    return df_result

  def _initSignificanceLevels(self):
    """
    Updates self.sl_min, self.sl_max: 
        self.avg, self.group, count, self.mutations
    """
    df_stat = self.df_stat.copy()
    def makeDF(val):
      """
      Constructs a 1 row dataframe selected by the value.
      """
      df = df_stat[df_stat[cn.AVG] == val].copy()
      df.reset_index(inplace=True)
      return df
    def assign(df, sl):
      """
      :param pd.DataFrame df: cn.COUNT, cn.GROUP, cn.AVG
      Assigns avg, count, group, mutations, group for SignificanceLevel
      """
      sl.avg = df.loc[0, cn.AVG]
      sl.group=  util.extendBin(
          bin(df.loc[0, cn.GROUP]), len(self._mutations))
      sl.count = df.loc[0, cn.COUNT]
      sl.mutations = self._mutations
    # Get mean values for each group
    maxval = df_stat[cn.AVG].max()
    df_max = makeDF(maxval)
    minval = df_stat[cn.AVG].min()
    df_min = makeDF(minval)
    # Assign the other values
    assign(df_min, self.sl_min)
    assign(df_max, self.sl_max)

  def calcMinMaxSLTStat(self):
    """
    Caculates the significance level for the minimum and maximum
    values of the means of groups obtained from the mutations
    using a T statistic.
    :param list-str mutations:
    :return float, float: significance level of min, max
    """
    num_obs = len(self._provider.df_y)
    mse = np.std(util.getFirstColumn(self._provider.df_y))
    mse = mse**2
    # Calculate the pooled standard deviation for each group
    df_stat = self.df_stat.copy()
    df_stat[cn.STD] = mse/df_stat[cn.COUNT]
    df_stat[cn.STD] = df_stat[cn.STD].apply(lambda v: np.sqrt(v))
    #
    def calcCDF(val, column):
      """
      Calculates the CDF for the val, putting the result in column.
      """
      df_stat[column] = val / df_stat[cn.STD]
      df_stat[column] = df_stat[column].apply(
        lambda v: t_dist.cdf(v, num_obs - 1))
    #
    calcCDF(self.sl_max.avg, CDF_PROB_MAX)
    calcCDF(self.sl_min.avg, CDF_PROB_MIN)
    # Calculate the significance level for the max
    sl_max = 1 - df_stat[CDF_PROB_MAX].prod()
    # Calculate the significance level for the min
    df_stat[CDF_PROB_MIN_ONE_MINUS] = df_stat[CDF_PROB_MIN].apply(
        lambda v: 1 - v)
    sl_min = 1 - df_stat[CDF_PROB_MIN_ONE_MINUS].prod()
    self.sl_max.sl_tstat = sl_max
    self.sl_min.sl_tstat = sl_min

  def calcMinMaxSLResample(self, num_replications=10000):
    """
    Caculates significance levels for the minimum and maximum
    values of the means of groups obtained from the mutations
    using resampling.
    :param list-str mutations:
    :param int num_replications:
    """
    df_stat = self.df_stat.copy()
    # Calculate the resample statistics
    group_sizes = util.getFirstColumn(
        self.df_group.groupby(cn.GROUP).count())
    y_vals = util.getFirstColumn(self._provider.df_y)
    df_resample = ud.generatePhenotypeData(y_vals,
        group_sizes, num_replications)
    df_mean = df_resample.groupby([cn.REPLICATION, cn.GROUP]).mean()
    df_mean.reset_index(inplace=True)
    #
    def findSL(val, is_max=True):
      """
      :return resampled significance level:
      """
      if is_max:
        df = df_mean.groupby(cn.REPLICATION).max()
      else:
        df = df_mean.groupby(cn.REPLICATION).min()
      vals = df[cn.VALUE].tolist()
      vals.sort()
      # Find position of val in vals
      if is_max:
        sl = (sum([0 if v < val else 1 for v in vals])*1.0)/len(vals)
      else:
        sl = (sum([0 if v > val else 1 for v in vals])*1.0)/len(vals)
      return sl
    #
    self.sl_min.sl_resample = findSL(self.sl_min.avg, is_max=False)
    self.sl_max.sl_resample = findSL(self.sl_max.avg, is_max=True)

  def getStatistics(self):
    """
    :return float, float: the non-nan values, else Bootstrap
    """
    def getStat(sl):
      if not np.isnan(sl.sl_resample):
        return sl.sl_resample
      else:
        return sl.sl_tstat
    #
    return getStat(self.sl_min), getStat(self.sl_max)
