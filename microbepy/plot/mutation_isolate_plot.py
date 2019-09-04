"""Provides plots of  mutations for Isolates and Lines."""

from microbepy.common import constants as cn
from microbepy.common.isolate import Isolate
from microbepy.common.study_context import StudyContext
from microbepy.data.model_data_provider import ModelDataProvider
from microbepy.plot.util_plot import PlotParms

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


COLORS = ['red', 'green', 'blue']
SPECIES = {cn.SPECIES_MIX_DVH: "DVH",
    cn.SPECIES_MIX_MMP: "MMP",
    None: "both"}
FONTSIZE_TITLE = 16
FONTSIZE_LABEL = 8


##################################################
# HELPER FUNCTIONS
##################################################
def _pruneColumns(df, species):
  columns = [c for c in df.columns if c[0] == species]
  for col in set(df.columns).difference(columns):
     del df[col]


##################################################
# Classes
##################################################
class MutationIsolatePlot(object):
  """
  Plot mutations with counts by isolates
  """

  def __init__(self, mutation_column=cn.GGENE_ID,
      provider=None, constraints=None,
      is_plot=True):
    """
    :param str mutation_column:
    :param ModelDataProvider provider:
    :param bool is_plot:
    """
    self._provider = provider
    self._mutation_column = mutation_column
    self._is_plot = is_plot
    if self._provider is None:
      context = StudyContext(depvar=cn.RATE,
          mutation_column=mutation_column)
      provider = ModelDataProvider(context, 
          constraints=constraints)
      provider.do()
    else:
      self._mutation_column = provider.context.mutation_column
    self.df_X = provider.df_X

  def _makeMutationDF(self, species):
    """
    :param str species:
    :return pd.DataFrame: columns are lines, index is mutation
    """
    index = self.df_X.index.tolist()
    df = self.df_X.copy()
    _pruneColumns(df, species)
    df[cn.LINE] = [Isolate.create(p[0]).line for p in index]
    df_result = pd.DataFrame()
    for idx, line in enumerate(set(df[cn.LINE])):
      df_result[line] = df[df[cn.LINE] == line].sum()
    df_result = df_result.drop([cn.LINE], axis=0)
    return df_result

  def plot(self, species):
    """
    Does a stacked bar plot of mutation frequency by isolate with 
    colors by line.
    :param str species:
    """
    df_plot = self._makeMutationDF(species)
    ax = df_plot.plot(kind='bar', stacked=True, figsize=(20,8), legend=None)
    ax.set_title("%s Mutations" % SPECIES[species], fontsize=FONTSIZE_TITLE)
    ax.set_ylabel("Isolate Count", fontsize=FONTSIZE_LABEL)
    ax.set_xlabel("Gene", fontsize=FONTSIZE_LABEL)
    plt.legend()
    if self._is_plot:
      plt.show()
