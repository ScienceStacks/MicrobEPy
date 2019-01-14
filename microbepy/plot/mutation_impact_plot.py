"""
Creates a plot showing the impact of a set of mutations.
"""

from microbepy.common import constants as cn
from microbepy.common.isolate import Isolate
from microbepy.common.mutation_context import MutationContext
from microbepy.common import util
from microbepy.data import util_data as ud
from microbepy.data.model_data_provider import ModelDataProvider
from microbepy.statistics.group_significance_level  \
    import GroupSignificanceLevel

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class MutationImpactPlot(object):

  def __init__(self, depvar, mutation_column,
      constraints=None, 
      provider=None,
      is_test=False):
    """
    :param str depvar: dependent variable
    :param str mutation_column:  Name of the column with mutations
    :param list-of-booleanFunction constraints:
    :param ModelDataProvider provider:
    """
    cls = self.__class__
    self.df_base = None
    #
    self._constraints = constraints
    self._is_test = is_test
    if provider is None:
      self._provider = ModelDataProvider(
          MutationContext(depvar, mutation_column), 
          constraints=constraints)
      self._provider.do()
    else:
      self._provider = provider
    self._depvar = self._provider.context.depvar
    self._mutation_column = self._provider.context.mutation_column
    self.df_X = self._provider.df_X
    self.df_y = self._provider.df_y

  def makeMutationStatisticDF(self, mutations):
    """
    Computes statistics for a group dataframe.
    :param list-str mutations:
    :return pd.DataFrame: cn.GROUP, cn.AVG, cn.STD, cn.COUNT
        The STD is the standard deviation of the mean
    """
    df = ud.makeMutationGroupDF(self.df_X, self.df_y, mutations)
    return genome_statistics.GenomeStatistics.makeGroupStatisticDF(
        df, cn.GROUP, cn.DEPVAR)
      
  def scatter(self, mutations, is_resample=False):
    """
    Construct a scatter plot for combinations of the
    mutations being present and absent.
    :param list-str mutations:
    :param bool is_resample: Use resample statistics
    """
    YLIM = [-3, 3]
    COLORS = 'bgry'
    group_sl = GroupSignificanceLevel(self._provider, mutations)
    if is_resample:
      group_sl.calcMinMaxSLResample()
    else:
      group_sl.calcMinMaxSLTStat()
    sl_min, sl_max = group_sl.getStatistics()
    # Limits on y values
    df = group_sl.df_group
    lines = df[cn.LINE].unique().tolist()
    if None in lines:
      lines.remove(None)
    color_dict = {l: COLORS[i] for i,l in enumerate(lines)}
    #
    for line in lines:
      plot_df = df[df[cn.LINE]==line]
      plt.scatter(plot_df[cn.GROUP], plot_df[cn.DEPVAR],
          c=color_dict[line])
    #
    df_stat = group_sl.makeGroupStatisticDF()
    df_stat = pd.DataFrame([r for _,r in df_stat.iterrows()
      if not np.isnan(r[cn.AVG])])
    plt.errorbar(df_stat[cn.GROUP], 
        df_stat[cn.AVG], 2*df_stat[cn.STD], marker='_', mfc='black',
             mec='black', ms=20, mew=4, ls='None', ecolor='black')
    # Add counts
    for group in df_stat[cn.GROUP]:
      df_val = df_stat[df_stat[cn.GROUP] == group]
      value = int(df_val[cn.COUNT].values[0])
      text = "(%d)" % value
      plt.text(group-0.2, YLIM[0]*0.9, text, fontsize=12)
    #
    plt.legend(lines)
    plt.ylim(YLIM)
    mutations = [str(m) for m in mutations]
    title = "Mutations: %s\nminSL: %1.4f, maxSL: %1.4f"   \
        % (', '.join(mutations), sl_min, sl_max)
    plt.title(title)
    plt.xlabel("Binary Encoded Mutations")
    plt.ylabel(self._depvar)
    if not self._is_test:
      plt.show()
    return df_stat
