"""Group by values of a dependent variable."""

import microbepy_init
import constants as cn
from dataframe_sorter import DataframeSorter
from model_data_provider import ModelDataProvider
from study_context import StudyContext
from range_constraint import RangeConstraintVector, RangeConstraint
import util

import copy
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from sklearn.cluster import KMeans


OUTPUT_DIRECTORY = util.getDataModelPath(None)
DEFAULT_RC = RangeConstraint(lower=-np.inf, upper=np.inf)
DEFAULT_RC_VECTOR = RangeConstraintVector({
    cn.RATE: DEFAULT_RC,
    cn.YIELD: DEFAULT_RC,
    })


class IsolateMutationStatistics(object):
  """
  Information about isolates and mutations.
  """

  def __init__(self, provider):
    self._provider = provider
    #
    self.mutations = None  # list-str
    self.isolate_pairs = None  # list-str
    self.ser_count = None  # pd.Series

  def do(self):
    self.mutations = self._provider.df_X.columns.tolist()
    self.isolate_pairs = self._provider.df_X.index.values
    self.ser_count = self._provider.df_X.sum()
    

class DepvarGroup(object):
  """
  Create and manage groups formed by being 
  in the same range of values of a dependent variable.
  Comparisons are made between isolates found within
  a constrained range and those that are unconstrained.
  """

  def __init__(self, 
      depvar, 
      mutation_column=cn.GGENE_ID,
      constraints=None,
      is_plot=True,
      output_directory=OUTPUT_DIRECTORY):
    """
    :param str depvar:
    :param str mutation_column:
    :param ModelDataProvider provider:
    :param list-BooleanFunction constraints:
    :param str output_directory:
    """
    self._depvar = depvar
    self._mutation_column = mutation_column
    self._constraints = constraints
    self._output_directory = output_directory
    self._is_plot = is_plot
    #
    self._provider = self.getData()  # Unconstrained by depvar
    self._statistics = IsolateMutationStatistics(self._provider)
    self._statistics.do()

  def getData(self, rc_vector=None):
    context = StudyContext(depvar=self._depvar,
        mutation_column=self._mutation_column)
    provider = ModelDataProvider(context,
        constraints=self._constraints,
        rc_vector=rc_vector)
    provider.do()
    return provider

  def makeStatisticsDF(self, lower, upper, steps):
    """
    Makes a statistics dataframe.
    :param float lower: lower end of range in units of std
    :param float upper: upper end of range in units of std
    :param float steps: number of steps between lower and upper
    :return pd.DataFrame:
        cn.LOWER: lower part of range of dependent variable
        cn.UPPER: upper part of range of dependent variable
        cn.CENTER: center of the range
        cn.ISOLATES: list of isolates in range
        cn.MUTATIONS: list of mutations in the range
        cn.FRACTION: Series with fraction of isolates with each mutation
    """
    step_size = (upper - lower)/(1.0*steps)
    result_dict = {
        cn.LOWER: [],
        cn.UPPER: [],
        cn.CENTER: [],
        cn.ISOLATES: [],
        cn.MUTATIONS: [],
        cn.FRACTION: [],
        }
    for step in range(steps):
      rc_vector = copy.deepcopy(DEFAULT_RC_VECTOR)
      new_lower = (lower + step_size*step)
      new_upper = new_lower + step_size
      rc_vector.update(self._depvar, RangeConstraint(
          lower=new_lower,
          upper=new_upper,
          ))
      provider = self.getData(rc_vector=rc_vector)
      statistics = IsolateMutationStatistics(provider)
      statistics.do()
      #
      result_dict[cn.LOWER].append(new_lower)
      result_dict[cn.UPPER].append(new_upper)
      result_dict[cn.CENTER].append((new_upper+new_lower)/2.0)
      result_dict[cn.ISOLATES].append(statistics.isolate_pairs)
      result_dict[cn.MUTATIONS].append(statistics.mutations)
      ser = statistics.ser_count / (1.0*len(statistics.isolate_pairs))
      result_dict[cn.FRACTION].append(ser)
    df_result = pd.DataFrame(result_dict)
    return df_result

  def makeHeatmap(self, steps, lower=-3, upper=3):
    """
    Creates a heatmap with depvar as x-axis and mutations as y.
    :param int steps:
    :param float lower: lower value for depvar
    :param float upper: upper value for depvar
    """
    max_major_labels = 5
    df_statistics = self.makeStatisticsDF(lower, upper, steps)
    # Construct a dataframe of fractions for each value of depvar
    df_plot = pd.DataFrame([r[cn.FRACTION] for _, r in 
        df_statistics.iterrows()])
    xvalues = df_statistics[cn.LOWER].tolist()
    df_plot.index = xvalues
    df_plot = df_plot.transpose()
    sorter = DataframeSorter(df_plot)
    df_plot = sorter.orderRows()
    # Plot
    fig = plt.figure(figsize=(16, 14))
    for idx, species in  \
        enumerate([cn.SPECIES_MIX_DVH, cn.SPECIES_MIX_MMP]):
      df_plot_species = df_plot.copy()
      ax = fig.add_subplot(1, 2, idx+1)
      mutations = [i for i in df_plot_species.index.values
          if i[0] == species]
      df_plot_species = df_plot_species.loc[mutations]
      plot = ax.pcolor(df_plot_species, cmap='jet')
      ax.set_yticks(np.arange(0.5, len(mutations), 1))
      ax.set_yticklabels(mutations, fontdict={'fontsize': 7})
      #
      num_minor_labels = max(1, len(xvalues)/max_major_labels)
      xlabels = [str(x) if i % num_minor_labels == 0 else "" 
           for i, x in enumerate(xvalues)]
      ax.set_xticks(np.arange(0.5, len(xvalues), 1))
      ax.set_xticklabels(xlabels)
      ax.set_xlabel(self._depvar)
      ax.set_title(cn.SPECIES_DICT[species])
      if idx == 1:
        fig.colorbar(plot, cmap='jet')
      else:
        ax.set_ylabel("Mutations")
    if self._is_plot:
      plt.show()
