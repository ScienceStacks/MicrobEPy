"""Module for common regression plots."""

import __init__
import constants as cn
from isolate import Isolate
import isolate_regression as ir
import isolate_model as im
from util_plot import PlotParms

import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np


################################################
# Classes
################################################
class RegressionPlot(object):
  #Regression plots of various sorts.

  def __init__(self, df, is_test=False):
    self._df = df
    self._is_test = is_test

  def estimateObserved(self, parms=None):
    """
    Plots estimated vs. observed values.
    :param pd.DataFrame df: columns: cn.ESTIMATE, cn.OBSERVED
    :param PlotParms params: plot parameters
    """
    if parms is None:
      parms = PlotParms()
    df_plot = self._df.copy()
    ax = df_plot.plot(kind='scatter', x=cn.OBSERVED, y=cn.ESTIMATE)
    ax.plot(self._df[cn.OBSERVED], self._df[cn.OBSERVED], 
        color='red')
    rsq = np.var(df_plot[cn.ESTIMATE])/np.var(df_plot[cn.OBSERVED])
    parms[cn.PLT_TITLE] = "%s: RSQ=%1.2f" % (parms[cn.PLT_TITLE], rsq)
    parms.do(is_plot=not self._is_test)

  def _calcMinMax(self):
    """
    Calculates the minimum and max values of the range for cn.ESTIMATE, cn.OBSERVED.
    :return (float, float):
    """
    lower = min(self._df[cn.ESTIMATE])
    lower1 = min(self._df[cn.OBSERVED])
    lower = min(lower, lower1)
    upper = max(self._df[cn.ESTIMATE])
    upper1 = max(self._df[cn.OBSERVED])
    upper = max(upper, upper1)
    return (lower, upper)

  def _makeLegendLine(self, cls_ir, isolate_pair, df_avg=None, 
      **kwargs):
    """
    :param IsolateRegression cls_ir:
    :param tuple-of-str isolate_pair:
    :param DataFrame df_avg:
    """
    epd_community = Isolate.create(isolate_pair[0]).getClonePairingID()
    legend_line = "%s" % epd_community
    if cls_ir == ir.AncestralPairingIsolateRegression:
      row = [r for _,r in df_avg.iterrows()
             if (r[cn.KEY_ISOLATE_DVH], r[cn.KEY_ISOLATE_MMP]) == isolate_pair][0]
      a_value = row[ir.COEF_DVH]
      b_value = row[ir.COEF_MMP]
      c_value = row[im.GAMMA]
      rsq_value = row[cn.RSQ]
      legend_line = "%s a=%2.4f, b=%2.4f, c=%2.4f, rsq=%1.2f" % (
          legend_line, a_value, b_value, c_value, rsq_value)
    else:
      legend_line = ""
    return legend_line

  def plotFilteredEstimatedObserved(self, cls_ir, depvar, max_std=3.0, 
      is_legend=True, leave_out_isolates=None, parms=None, **kwargs):
    """
    Plots estimated vs. observed values that are filtered.
    :param IsolateRegression-Subclass cls_ir:
    :parm str depvar:
    :param float max_std: maximum standard deviation for the residuals
    :param list-of-str-tuple leave_out_isolates:
    :param dict kwargs: other parameters needed by cls_ir
    """
    if parms is None:
      parms = PlotParms()
      parms[cn.PLT_FIGSIZE] = (12, 8)
    def keepIsolates(row):
      if leave_out_isolates is None:
        return True
      if (row[cn.KEY_ISOLATE_DVH], row[cn.KEY_ISOLATE_MMP]) in leave_out_isolates:
        return False
      return True
    #   
    # Filter based on standard deviation of residuals
    cultures = cls_ir.getSmallResidualCultures(max_std)
    isPermittedRow = lambda r: keepIsolates(r) and (r[cn.KEY_CULTURE] in cultures)
    dfs = cls_ir.makeEstimateDFS(depvar=depvar, isPermittedRow=isPermittedRow, **kwargs)
    self._df = dfs[cn.RESIDUAL].copy()
    self._df.reset_index(inplace=True)
    cultures = self._df[cn.KEY_CULTURE].tolist()
    # Construct the title
    df_estimate = dfs[cn.AVG].copy()
    # Plot
    groups = self._df.groupby([cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP]).groups
    isolate_pairs = list(groups.keys())
    colors = cm.rainbow(np.linspace(0, 1, len(isolate_pairs)))
    legend_labels = []
    for isolate_pair, color in zip(isolate_pairs, colors):
      indices = groups[isolate_pair]
      xvalues = self._df.loc[indices, cn.OBSERVED].tolist()
      yvalues = self._df.loc[indices, cn.ESTIMATE].tolist()
      label = self._makeLegendLine(cls_ir, 
          isolate_pair, df_avg=df_estimate)
      plt.scatter(xvalues, yvalues, color=color, label=label)
    point = self._calcMinMax()
    plt.plot(point, point, 'r-', label="unit slope")
    #
    parms[cn.PLT_TITLE] = "%s %s: RSQ=%1.3f, max_std=%2.2f" % (
        parms[cn.PLT_TITLE], depvar, dfs[cn.RSQ], max_std)
    if is_legend:
      parms[cn.PLT_LEGEND] = None
    parms.do(is_plot=not self._is_test)
