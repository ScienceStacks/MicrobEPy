"""
   Quantifies the difference in occurrence of one
   (or more mutations) between two sets of
   isolates.

   TODO:
    1. Use DualDataProvider and implement getPoints locally so
    have same data plotted as is analyzed
"""


from microbepy.common import constants as cn
from microbepy.common import util
from microbepy.common.range_constraint import RangeConstraint
from microbepy.data import util_data as ud
from microbepy.data import model_data_provider
from microbepy.plot.phenotype_plot import PhenotypePlot
from microbepy.plot.util_plot import PlotParms

from collections import namedtuple
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


LOW = 0.0
# Depvar is less than 0
RC_LOW = RangeConstraint(lower=-10, upper=-LOW)  # Units of standard deviation
# Depvar is greater than 0
RC_HIGH = RangeConstraint(lower=LOW, upper=10)


######################################
# CLASSES
######################################

# Knobs are the protrusions from points
Knob = namedtuple('Knob', ['angle', 'color'])

Point = namedtuple('Point', ['x', 'y'])


class MutationDifferential(object):

  def __init__(self, depvar, mutation_column, 
      rc_low=RC_LOW,
      rc_high=RC_HIGH,
      is_median=False,
      **kwargs
      ):
    """
    The sets of isolates are specified by a range of values
    of the dependent variable.
    :param str depvar:
    :param str mutation_column:
    :param RangeConstraint rc_low: Defines the low isolates
    :param RangeConstraint rc_high: Defines the high isolates
    :param Boolean is_median: use the median as 0 if True
    :param dict **kwargs: arguments to provider
    """
    def shiftRC(rc, shift):
      return RangeConstraint(rc.lower + shift, rc.upper + shift)
    #
    self._depvar = depvar
    self._mutation_column = mutation_column
    self._rc_low = rc_low
    self._rc_high = rc_high
    if not "constraints" in kwargs:
      kwargs["constraints"] = None
    self._constraints = kwargs["constraints"]
    # Obtain the data
    self.provider = model_data_provider.makeTransformedData(
        data_cls=model_data_provider.ModelDataDualProvider,
        mutation_column=self._mutation_column,
        **kwargs)
    self.provider.do()
    self.df_y = self.provider.df_ys[depvar].rename(
        columns={depvar: cn.VALUE})
    self.df_X = self.provider.df_X
    # Adjust the range constraints if the median is used
    if is_median:
      median = np.median(self.df_y[cn.VALUE])
      self._rc_low = shiftRC(self._rc_low, median)
      self._rc_high = shiftRC(self._rc_high, median)
    # Mutation counts for low range isolates
    self._ser_low, self._tot_low = self._makeCounts(self._rc_low)
    # Mutation counts for high range isolates
    self._ser_high, self._tot_high = self._makeCounts(self._rc_high)

  def _makeCounts(self, range_constraint):
    """
    :param RangeConstraint range_constraint:
    :return pd.Series, int: 
        pd.Series: index is mutation; values is count
        int: count of isolates
    """
    sel = range_constraint.isSatisfiedRows(self.df_y)
    df = self.df_X.loc[sel]
    return df.sum(axis=0), len(sel)

  @staticmethod
  def _calcProb(num_low, num_high, tot_low, tot_high):
    """
    Caculates the probability of having a mutation
    that occurs num_low times in the low isolates and
    num_high times in the high isolates.
    :param int num_low: Number in low slots
    :param int num_high: Number in high slots
    :param int tot_low: total low slots
    :param int tot_high: total high slots
    :return float:
    """
    tot = tot_low + tot_high
    result = util.nCr(tot_low, num_low)
    result *= util.nCr(tot_high, num_high)
    result /= 1.0*util.nCr(tot, num_low + num_high)
    return result

  def _calcSL(self, num_low, num_high):
    """
    Calculates the significance level for a mutation
    that has num_low mutations in the low isolates in
    num_high in the high isolates. Two cases are
    considered: (1) The mutation is unusual in its
    occurrence in the high isolates; (2) the mutation
    is unusual in its absence in the high isolates.
    :param int num_low:
    :param int num_high:
    :return float:
    """
    # Significance for mutation present in high
    #   Find the tail by moving lows to high
    prob_tail_high = 0
    upper = num_low + 1
    for num in range(0, upper):
      n_high = num_high + num
      if n_high > self._tot_high:
        break
      n_low = num_low - num
      prob_tail_high += MutationDifferential._calcProb(
          n_low, n_high, self._tot_low, self._tot_high)
    # Find the tail by moving highs
    prob_tail_low = 0
    upper = num_high + 1
    for num in range(0, upper):
      n_high = num_high - num
      n_low = num_low + num
      if n_low > self._tot_low:
        break
      prob_tail_low += MutationDifferential._calcProb(
          n_low, n_high, self._tot_low, self._tot_high)
    #
    sl = min(prob_tail_high, prob_tail_low)
    return sl

  def makeDF(self):
    """
    A significance level is calculated for preferentially
    occurring the in the high (low) set.
    :param int num_replications:
    :return pd.DataFrame: sorted by increasing value of cn.SIGLVL
        indexed by cn.MUTATION
        cn.SIGLVL - float significance level 
            for the occurrence of the mutation being preferentially in the high (low)
            group.
        cn.VALUE - joint significance level of mutations with
                   up this mutation in list
        cn.COUNT1 - number low
        cn.COUNT2 - number high
    """
    # Calculate significance levels
    sls = [self._calcSL(int(l), int(h)) for l, h 
        in zip(self._ser_low, self._ser_high)]
    # Construct the dataframe
    df_result = pd.DataFrame({
        cn.MUTATION: self.df_X.columns.tolist(),
        cn.COUNT1: self._ser_low,
        cn.COUNT2: self._ser_high,
        cn.SIGLVL: sls,
        })
    df_result = df_result.set_index(cn.MUTATION)
    df_result = df_result.sort_values(cn.SIGLVL)
    # Compute joint significance level
    joint_prob = 1
    values = []
    for _, row in df_result.iterrows():
      joint_prob *= (1 - row[cn.SIGLVL])
      values.append(1 - joint_prob)
    df_result[cn.VALUE] = values
    #
    return df_result

  def scatterKnob(self, is_plot=True, parms=PlotParms(), max_sl=0.1):
    """
    Plots isolates in phenotype space with knob marks that indicate
    mutations that are statistically significant.
    :param boolean is_plot: plot if True
    :param PlotParms parms: plot parameters
    :param float max_sl:
    
    TODO:
      1. Partition graph with dotted line along the cutting axis
      2. Legend with mutations
    """
    def plotKnob(point, knob, length, y_adjust):
     """
     Plots the knob for the point
     :param Point point:
     :param Knob knob:
     :param float length: length of knob
     :param float y_adjust: adjust to y because of x, y scales
     """
     # find the end point
     endy = y_adjust*length * math.sin(
         math.radians(knob.angle)) + point.y
     endx = length * math.cos(
         math.radians(knob.angle)) + point.x
     plt.plot([point.x, endx], [point.y, endy], 
         color=knob.color, linewidth=3)
    #
    knobs = [
        Knob(angle=0, color='r'),
        Knob(angle=90, color='g'),
        Knob(angle=180, color='b'),
        Knob(angle=270, color='y'),
        Knob(angle=45, color='m'),
        Knob(angle=135, color='k'),
        ]
    num_knobs = len(knobs)
    # Construct the base plot
    plot = PhenotypePlot(provider=self.provider,
         constraints=self._constraints, is_plot=False)
    parms.setattr(cn.PLT_XLIM, [-2, 2])
    parms.setattr(cn.PLT_YLIM, [-3.5, 3.5])
    plot.scatter(colors='kkk', legend=[], parms=parms)
    df_point = plot.getPoints()
    # Add knobs for significant mutations
    df_mutation = self.makeDF()  # Significant mutations
    df_sig = df_mutation[df_mutation[cn.SIGLVL] < max_sl].copy()
    mutations = df_sig.index.tolist()[:num_knobs]
    siglvls = ["%1.4f" % s for s in df_sig[cn.SIGLVL]]
    # 
    length = 0.08
    y_adjust = parms[cn.PLT_YLIM][1]/(1.0*parms[cn.PLT_XLIM][1])
    for idx, mutation in enumerate(mutations):
      ser = self.df_X[mutation]
      pairs = [i for v,i in zip(ser.tolist(), ser.index.tolist())
          if v == 1]
      df_isolates = pd.DataFrame(
          [r for p,r in df_point.iterrows() if p in pairs]
          )
      points = [Point(x=r[cn.RATE], y=r[cn.YIELD])
          for _, r in df_isolates.iterrows()]
      for point in points:
        plotKnob(point, knobs[idx], length, y_adjust)

    # Add the line dividing the isolates being compared
    vals = [self._rc_low.upper, self._rc_high.lower]
    if self._depvar == cn.RATE:
      xs = vals
      ys = parms[cn.PLT_YLIM]
    else:  # cn.YIELD
      ys = vals
      xs = parms[cn.PLT_XLIM]
    plt.plot(xs, ys, dashes=[10, 5, 20, 5], linewidth=1, 
        color='black')
    # Add the legend
    if not cn.PLT_LEGEND in parms.keys():
      entries = []
      for idx, mutation in enumerate(mutations):
        label = "%s (%s)" % (mutations[idx], siglvls[idx])
        entries.append(mpatches.Patch(color=knobs[idx].color, label=label))
    else:
      entries = parms[cn.PLT_LEGEND]
    # Position the legend to get it out of view
    plt.legend(handles=entries, loc='lower left')
    # Show the plot
    parms.do(is_plot=is_plot)
