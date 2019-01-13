"""
A collection of mutation groups, where each group consists of
mutations that occur in the same isolates.
"""

from microbepy.common import constants as cn
from microbepy.common.dataframe_sorter import DataframeSorter
from microbepy.common import group_collection
from microbepy.data import model_data_provider
from microbepy.plot.util_plot import PlotParms
from microbepy.common import util

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


##################### CLASSES ##########################
class MutationCollection(group_collection.GroupCollection):

  def __init__(self, initial_groups=None, prefix="", 
      is_set_group_label=True, is_plot=True):
    """
    If no mutation_source is specified, then calculates a default
    MutationDifferential.
    :param list-Group initial_groups:
    :param bool is_plot: plots if True
    """
    initial_groups = util.setNoneList(initial_groups)
    super(self.__class__, self).__init__(initial_groups=initial_groups)
    if is_set_group_label:
      self.setGroupLabels(prefix=prefix)
    self._is_plot = is_plot

  @staticmethod
  def calcSL(num_isolates, num_mutations, num_occurrences):
    """
    :param int num_isolates: Number of isolates in line
    :param int num_mutations: Mutations that co-occur
    :param int num_occurrences: Number of isolates in which
        the mutations co-occur
    :return float: probability
    """
    num_permute = util.nCr(num_isolates, num_occurrences)
    denom = num_permute**num_mutations
    return 1.0/denom

  @classmethod
  def makeGroups(cls, provider, species=cn.SPECIES_MIX_DVH):
    """
    :param ModelDataProvider provider:
    :return list-Group
    :param str species:
    """
    groups = []
    num_isolates = len(provider.df_X)
    for column in provider.df_X.columns:
      if group_collection.ELEMENT_SEPARATOR in column:
        mutations = column.split(group_collection.ELEMENT_SEPARATOR)
        selected_mutations = [m for m in mutations if m[0] == species]
        if len(selected_mutations) < 2:
          continue
        group = group_collection.Group(selected_mutations)
        num_occurrences = provider.df_X[column].sum()
        prob = cls.calcSL(num_isolates, group.len(), num_occurrences)
        group.value = -np.log10(prob)  # Number 0s in sign. level
        groups.append(group)
    return groups

  @classmethod
  def makeMutationCollectionForLine(cls, line=None, 
      mutation_column=cn.GGENE_ID, species=cn.SPECIES_MIX_DVH,
      **kwargs):
    """
    Calculates the GroupCollection for a line. Computes
    the -log10(significance level) of each group.
    :param str line: if None, do all lines
    :param str species:
    :param dict kwargs: Arguments for MutationCollection constructor
    :return MutationCollection:
    """
    def processLine(line):
      constraints = [ lambda r: r[cn.LINE] == line]
      provider = model_data_provider.makeTransformedData(
          constraints=constraints, mutation_column=mutation_column)
      groups = cls.makeGroups(provider, species=species)
      collection = MutationCollection(groups, prefix=line,
          is_set_group_label=True, **kwargs)
      return collection
    #
    if line is not None:
      return processLine(line)
    #
    collection = MutationCollection([])  # Empty collection
    for line in cn.LINE_CIS:
      new_collection = processLine(line)
      collection = collection.unionDisjoint(new_collection,
          is_set_group_label=False, **kwargs)
    return collection

  def plot(self, values_dict={}, parms=PlotParms()):
    """
    Plots groups on x-axis and group objects on y-axis.
    Clusters like groups and like objects.
    :param PlotParms parms:
    
    TODO
      1. Have labels for groups so can distinguish lines?
         Ex: HA2-1, UE3-4
      2. Groups may have a significance level? Use -log(SL),
         which is the number of 0's following the decimal point.
         Compute significance levels for each group.
    """
    fig = plt.figure(figsize=parms[cn.PLT_FIGSIZE])
    parms[cn.PLT_XLABEL] = "Mutation Group"
    parms[cn.PLT_YLABEL] = "Mutation"
    ax = fig.add_subplot(1, 1, 1)
    sorter = DataframeSorter(self.makeValueDF())
    df_plot = sorter.orderRows()
    # Order the columns based on similarity in mutations present
    plot = ax.pcolor(df_plot, cmap='jet')
    labels = df_plot.columns.tolist()
    ax.set_xticks(np.arange(0.5, len(labels)))
    ax.set_xticklabels(labels, rotation=90)
    mutations = df_plot.index.tolist()
    ax.set_yticks(np.arange(0.5, len(mutations)))
    ax.set_yticklabels(mutations)
    # fig.colorbar(plot, cmap='jet', boundaries=range(20))
    fig.colorbar(plot, cmap='jet')
    parms.do(is_plot=self._is_plot)
