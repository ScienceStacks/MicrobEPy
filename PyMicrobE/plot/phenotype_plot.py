"""
Plots values in phenotype space: cn.RATE, cn.YIELD
"""

import __init__
import constants as cn
import util
import util_data as ud
from model_data_provider import ModelDataDualProvider, ModelDataProvider
from group_significance_level import GroupSignificanceLevel
from study_context import StudyContext
from util_plot import PlotParms

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class PhenotypePlot(object):

  def __init__(self, constraints=None, 
      provider=None, is_plot=True):
    """
    :param list-of-booleanFunction constraints:
    :param ModelDataDualProvider provider:
    """
    cls = self.__class__
    #
    self.provider = provider
    if provider is None:
      # Arbitrary choice of mutation_column since only want df_y
      self.provider = ModelDataDualProvider(cn.GGENE_ID,
          constraints=constraints)
      self.provider.do()
    self.lines = ModelDataProvider.getLinesForRows(self.provider.df_X)
    self._is_plot = is_plot
      
  def scatter(self, parms=PlotParms(), is_errorbars=False,
      colors='bgry', legend=None):
    """
    Construct a scatter plot of RATE vs. YIELD
    coloring by line.
    :param PlotParms parms: plots specifics
    :param bool is_errorbars: include error bars
    :param str colors: to plot
    :param list-str legend:
    """
    XLIM = [-2, 2]
    YLIM = [-3.5, 3.5]
    unique_lines = list(set(self.lines))
    color_dict = {l: colors[i] for i,l in enumerate(unique_lines)}
    #
    for line in unique_lines:
      sel = [line == l for l in self.lines]
      xs = self.provider.df_ys[cn.RATE].loc[sel, cn.VALUE]
      ys = self.provider.df_ys[cn.YIELD].loc[sel, cn.VALUE]
      plt.scatter(xs, ys, c=color_dict[line])
      # Error bars
      if is_errorbars:
        yerr = [v for v in 
            self.provider.df_y_stds[cn.YIELD].loc[sel, cn.VALUE]]
        xerr = [v for v in 
            self.provider.df_y_stds[cn.RATE].loc[sel, cn.VALUE]]
        plt.errorbar(xs, ys, xerr=xerr, yerr=yerr, c=color_dict[line],
            fmt='o', elinewidth=0.75)
    #
    if legend is None:
      if len(unique_lines) > 1:
        plt.legend(unique_lines)
      else:
        parms[PLT_TITLE] = unique_lines[0]
    else:
      plt.legend(legend)
    parms[cn.PLT_XLABEL] = cn.RATE
    parms[cn.PLT_YLABEL] = cn.YIELD
    parms.do(is_plot=self._is_plot)

  def getPoints(self):
    """
    :return pd.DataFrame: cn.RATE, cn.YIELD indexed by ISOLATE_PAIR
    """
    df = pd.DataFrame({
      cn.RATE: self.provider.df_ys[cn.RATE][cn.VALUE].tolist(),
      cn.YIELD: self.provider.df_ys[cn.YIELD][cn.VALUE].tolist(),
      })
    df.index = self.provider.df_ys[cn.RATE].index
    return df
