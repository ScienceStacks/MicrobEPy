"""Provides plots of  mutations for Isolates and Lines."""

from microbepy.common import constants as cn
from microbepy.common.dataframe_sorter import DataframeSorter
from microbepy.common.isolate import Isolate
from microbepy.common.study_context import StudyContext
from microbepy.common import util
from microbepy.correlation import genome_correlation
from microbepy.data.model_data_provider import ModelDataProvider
from microbepy.plot.util_plot import PlotParms

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COLORS = ['red', 'green', 'blue']
SPECIES = {cn.SPECIES_MIX_DVH: "DVH",
    cn.SPECIES_MIX_MMP: "MMP"}
FONTSIZE_TITLE = 16
FONTSIZE_LABEL = 14
MAX_LINES = 13
MIN_FRACTION = 0.25
MAX_SIGLVL = 0.01
COLORBAR_MIN = 1.0
COLORBAR_MAX = 4.0


########################################################################
class MutationPlot(object):
  """
  Common code
  """

  @staticmethod
  def _pruneColumns(df, species):
    columns = [c for c in df.columns if c[0] == species]
    for col in set(df.columns).difference(columns):
       del df[col]


########################################################################
class MutationIsolatePlot(MutationPlot):
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
    self.__class__._pruneColumns(df, species)
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


#####################################################################
# TODO: Heat plot for significance levels, sorted by similarity
class MutationLinePlot(MutationPlot):
  """
  Plot mutations by occurrences within Lines.
  """

  def __init__(self, mutation_column=cn.GGENE_ID,
      is_plot=True):
    """
    :param str mutation_column:
    :param bool is_plot:
    """
    self._mutation_column = mutation_column
    self._is_plot = is_plot

  def getTransfers(self):
    """
    :return list-int: transfers sorted low to high
    """
    query = '''
        select distinct transfer from genotype
        where key_isolate like '%%*.*.*.*' 
            and is_an_mutation = 0
        '''
    df = util.readSQL(query)
    result = [int(s) for s in df[cn.TRANSFER].unique()]
    result.sort()
    return result

  def getLines(self, species=None):
    """
    Obtains the lines present for one or both species.
    :param str species:
    :return list-int: transfers sorted low to high
    """
    def query(species):
      query = '''
          select distinct line from genotype
          where is_an_mutation = 0
              and key_mutation like '%s%%'
          ''' % species
      df = util.readSQL(query)
      return df[cn.LINE].unique()
    #
    if species is not None:
      result = query(species)
    else:
      result_dvh = query(cn.SPECIES_MIX_DVH)
      result_mmp = query(cn.SPECIES_MIX_MMP)
      result = list(set(result_dvh).intersection(result_mmp))
    result.sort()
    return result

  def _makeLineDF(self, species=None, permitted_mutations=None,
       transfer=cn.TRANSFER_DEFAULT, **kwargs):
    """
    :param str species:
    :param list-str permitted_mutations:
    :params dict kwargs: parameters passed to select mutations
    :return pd.DataFrame: columns are mutation with values of freq
        indexed by line
    """
    # Get the data
    query = '''
    select distinct line, key_isolate, ggene_id, freq 
      from genotype
      where 
          key_isolate like '%%*.*.*.*' 
          and is_an_mutation = 0 
          and transfer = %s
      order by line, transfer, key_isolate, ggene_id
    ''' % str(transfer)
    df = util.readSQL(query)
    sel = [Isolate.create(i).epd_id == cn.ISOLATE_DEFAULT
        for i in df[cn.KEY_ISOLATE]]
    df = df[sel].copy()
    del df[cn.KEY_ISOLATE]
    # Construct the matrix with all relevant lines
    df_matrix = util.makeMatrix(df, row_name=self._mutation_column,
        column_name=cn.LINE, value_name=cn.FREQ)
    lines = self.getLines()
    for line in set(lines).difference(df_matrix.columns):
      df_matrix[line] = 0.0
    lines.sort()
    df_matrix = df_matrix[lines]
    # Select the species
    if species is not None:
      df_transpose = df_matrix.transpose()
      self.__class__._pruneColumns(df_transpose, species)
      df_matrix = df_transpose.transpose()
    # Adjust for permitted mutations
    if permitted_mutations is None:
      permitted_mutations = set(self._getFrequentMutations(
          species=species, **kwargs))
    for mutation in df_matrix.index:
      if not mutation in permitted_mutations:
        df_matrix = df_matrix.drop(mutation)
    for mutation in permitted_mutations:
      if not mutation in df_matrix.index:
        row = {k: 0 for k in df_matrix.columns}
        df_matrix.loc[mutation] = row
    df_matrix = df_matrix.loc[permitted_mutations]
    # Convert from percent to fractions
    df_matrix = df_matrix.applymap(lambda v: 0.01*v)
    #
    return df_matrix

  def _orderMutations(self, species=None, **kwargs):
    """
    Orders mutations based on similarities in occurrences in lines.
    :param str species:
    :param dict kwargs: optional arguments
    :return list-str:
    """
    transfers = self.getTransfers()
    df_lines = pd.DataFrame()
    for transfer in transfers:
      df = self._makeLineDF(species=species, permitted_mutations=None,
          transfer=transfer, **kwargs)
      for column in df.columns:
        new_column = "%s_%d" % (column, transfer)
        df_lines[new_column] = df[column]
    #
    sorter = DataframeSorter(df_lines)
    df_sort = sorter.orderRows()
    return df_sort.index.tolist()

  def plotTransfers(self, species=None, 
        parms=PlotParms(is_initialize=False), 
        is_cluster_mutations=True,
        **kwargs):
    """
    Does a stacked bar plot of mutation frequency for all transfers.
    :params str species:
    :param bool is_cluster_mutations: Group similar mutations together
    :return pd.DataFrame: row=mutation, col=line + transfer, value is fraction
    """
    if is_cluster_mutations:
      permitted_mutations = self._orderMutations(
          species=species, **kwargs)
    else:
      permitted_mutations = self._getFrequentMutations(
          species=species)
    transfers = self.getTransfers()
    num_transfers = len(transfers)
    fig, axes = plt.subplots(nrows=num_transfers, ncols=1)
    dfs = []
    for idx, transfer in enumerate(transfers):
      parms[cn.PLT_YTICKLABELS] = True
      if species is None:
        parms[cn.PLT_TITLE] = "%d" % transfer
      else:
        parms[cn.PLT_TITLE] = "%s, %d" % (species, transfer)
      if idx == 0:
        parms[cn.PLT_YLABEL] = True
      else:
        parms[cn.PLT_YLABEL] = False
      if idx < num_transfers - 1:
        parms[cn.PLT_LEGEND] = False
        parms[cn.PLT_XLABEL] = False
        parms[cn.PLT_XTICKLABELS] = False
      else:
        parms[cn.PLT_LEGEND] = True
        parms[cn.PLT_XLABEL] = True
        parms[cn.PLT_XTICKLABELS] = True
      df = self.plotLine(species, transfer, parms=parms, is_plot=False,
          ax=axes[idx], permitted_mutations=permitted_mutations,
          **kwargs)
      df[cn.TRANSFER] = transfer
      dfs.append(df)
    if self._is_plot:
      plt.show()
    return pd.concat(dfs)
    

  def plotLine(self, species, transfer, 
      parms=PlotParms(is_initialize=False),
      is_plot=None, ax=None, permitted_mutations=None, **kwargs):
    """
    Does a stacked bar plot of mutation frequency by line
    with colors
    :params str species: If None, do both
    :params int transfer:
    :params PlotParms parms:
    :params Axis ax: axis to use in plot
    :param list-str permitted_mutations: to use and how they
       are ordered if None, then use alphabetical order
    :params dict kwargs: parameters passed to select mutations
    :return pd.DataFrame: row=mutation, col=line, value is fraction
    """
    if is_plot is None:
      is_plot = self._is_plot
    parms.setTrueIfAbsent(cn.PLT_XLABEL)
    parms.setTrueIfAbsent(cn.PLT_XTICKLABELS)
    #
    df_plot = self._makeLineDF(species=species, 
        permitted_mutations=permitted_mutations,
        transfer=transfer, **kwargs)
    # Do the plot
    if not cn.PLT_FIGSIZE in parms:
      parms[cn.PLT_FIGSIZE] = (12, 8)
    if ax is None:
      ax = df_plot.plot(kind='bar', stacked=True, 
          figsize=parms[cn.PLT_FIGSIZE], legend=None)
    else:
      df_plot.plot(kind='bar', stacked=True, 
          legend=None, ax=ax, figsize=parms[cn.PLT_FIGSIZE])
    ax.set_xlabel("", fontsize=FONTSIZE_LABEL)  # Eliminate implicit label
    if parms.isFalse(cn.PLT_XTICKLABELS):
      labels = ax.get_xticklabels()
      new_labels = np.repeat("", len(labels))
      ax.set_xticklabels(new_labels)
    if parms.isFalse(cn.PLT_YTICKLABELS):
      labels = ax.get_yticklabels()
      new_labels = np.repeat("", len(labels))
      ax.set_yticklabels(new_labels)
    if cn.PLT_TITLE in parms:
      title = parms[cn.PLT_TITLE]
    else:
      title = "%s Mutations" % SPECIES[species]
    xpos = int(len(df_plot)*0.5)
    ypos = MAX_LINES - 3
    ax.text(xpos, ypos, title, fontsize=FONTSIZE_TITLE)
    #ax.set_title(title, fontsize=FONTSIZE_TITLE)
    ax.set_ylim([0, MAX_LINES])
    if parms.isTrue(cn.PLT_YLABEL):
      ax.set_ylabel("Fraction", fontsize=FONTSIZE_LABEL)
    if parms.isTrue(cn.PLT_XLABEL):
      ax.set_xlabel(self._mutation_column, fontsize=FONTSIZE_LABEL)
    if parms.isTrue(cn.PLT_LEGEND):
      ax.legend(loc=(1,2))
      #ax.legend()
    if is_plot:
      plt.show()
    return df_plot

  def _getFrequentMutations(self, species=cn.SPECIES_MIX_DVH, 
      min_lines=2):
    """
    :param str species:
    :param int min_lines: minimum number of lines in which
        mutation occurs
    :return list-str: 
    """
    query = '''
        select distinct %s, count(distinct line) as %s
          from genotype
          where is_an_mutation = 0
          group by %s
          order by count(distinct line) DESC
        ''' % (self._mutation_column, cn.COUNT, self._mutation_column)
    df = util.readSQL(query)
    mutations = [r[self._mutation_column] for _, r in df.iterrows()
        if r[cn.COUNT] >= min_lines]
    if species is not None:
      mutations = [m for m in mutations if m[0] == species]
    mutations.sort()
    return mutations

  def _makeMutationSiglvlMatrix(self, species=None,
       transfer=cn.TRANSFER_DEFAULT, 
       other_transfer=None, min_fraction=MIN_FRACTION,
       **kwargs):
    """
    Creates a significance level matrix for mutations.
    :param str species:
    :param int transfer: transfer time for row mutations
    :param int other_transfer: transfer time for column mutations
    :param float min_fraction: minimum fractional occurrence of
        a mutation within a line for it to be considered
    :params dict kwargs: parameters passed to select mutations
    :return pd.DataFrame: row index and columns are mutations
    """
    def makeDF(transfer):
      df_line = self._makeLineDF(species=species, transfer=transfer,
          **kwargs)
      df_binary = df_line.applymap(
          lambda v: 1.0 if v > min_fraction else 0)
      return df_binary.transpose()
    #
    if other_transfer is None:
      other_transfer = transfer
    #
    df_binary_rows = makeDF(transfer)
    df_binary_columns = makeDF(other_transfer)
    df_matrix = genome_correlation.makeSiglvlDF(df_binary_rows,
        df_other=df_binary_columns)
    return df_matrix

  def _plotSiglvlDF(self, transfer=cn.TRANSFER_DEFAULT,
      other_transfer=None,
      max_siglvl=MAX_SIGLVL,
      **kwargs):
    """
    Constructs a the dataframe used for heatmap.
    :param int transfer:
    :param float max_siglvl:
    :return pd.DataFrame: mutions, mutations,
        values are -log10 significance level
    """
    df_matrix = self._makeMutationSiglvlMatrix(transfer=transfer,
        other_transfer=other_transfer, **kwargs)
    sorter = DataframeSorter(df_matrix)
    df_sort = sorter.orderBoth()
    #
    df_transformed = df_sort.applymap(lambda v: np.log10(v))
    df_transformed = df_transformed.applymap(lambda v: -v)
    ubound = -np.log10(max_siglvl)
    df_plot = df_transformed.applymap(
        lambda v: np.nan if v < ubound else v)
    sorter = DataframeSorter(df_plot)
    df_plot = sorter.deleteNanRowsAndColumns()
    return df_plot

  def plotSiglvls(self, is_time_lag=False, 
      parms=PlotParms(), **kwargs):
    """
    Does a subplots of mutation correlation significance levels.
    :param bool is_time_lag: construct time lag subplots
    :param dict kwargs: non-transfer parameters passed to next level
    :return dict: key is pair of transfers, value is data_frame
    """
    NCOLS = 3
    NPLOTS = 9
    plot_pos = {1:1, 2:3, 3:4, 4:6, 5: 7}
    transfers = self.getTransfers()
    if is_time_lag:
      pairs = [p for p in zip(transfers[0:-1], transfers[1:])]
    else:
      pairs = [p for p in zip(transfers[:-1], transfers[:-1])]
    #
    nrows = 2 if (len(pairs) == 4) else 3
    fig = plt.figure(figsize=parms[cn.PLT_FIGSIZE])
    result = {}
    for idx, pair in enumerate(pairs):
      idx += 1
      ax = fig.add_subplot(nrows, NCOLS, plot_pos[idx])
      if idx < len(pairs):
        is_plot = False
      else:
        is_plot = True
      if idx in [1, 2, 5]:
        parms[cn.PLT_XAXISTICKTOP] = True
      else:
        parms[cn.PLT_XAXISTICKTOP] = False
      if idx == 4:
        parms[cn.PLT_COLORBAR] = True
      else:
        parms[cn.PLT_COLORBAR] = False
      df = self.plotSiglvl(transfer=pair[0], other_transfer=pair[1],
          fig=fig, ax=ax, parms=parms, is_plot=is_plot, **kwargs)
      result[pair] = df
    return result

  def plotSiglvl(self, transfer=cn.TRANSFER_DEFAULT,
      other_transfer=None,
      max_siglvl=MAX_SIGLVL,
      ax=None,
      fig=None,
      parms=PlotParms(),
      is_plot=None):
    """
    Constructs a heatmap of the mutation coocurrence significance
    levels.
    :param int transfer:
    :param int other_transfer: Allow comparisons across time
    :param Matplotlib.Axes ax:
    :param PlotParms parms: Parameters for the plot
    :param bool is_plot: Overrides constructor plotting directive
    :return pd.DataFrame: columns, rows are mutations
    """
    def makeLabel(transfer, column):
      return "%d-%s" % (transfer, column)
    #
    if is_plot is None:
      is_plot = self._is_plot
    elif not self._is_plot:
      is_plot = self._is_plot
    #
    df_plot = self._plotSiglvlDF(transfer=transfer,
        other_transfer=other_transfer,
        max_siglvl=max_siglvl)
    mutations = df_plot.columns.tolist()
    # Do the plot
    if not cn.PLT_COLORBAR in parms:
      parms[cn.PLT_COLORBAR] = True
    if other_transfer is None:
        other_transfer = transfer
    if ax is None:
      if fig is None:
        fig = plt.figure(figsize=parms[cn.PLT_FIGSIZE])
      ax = fig.add_subplot(1, 1, 1)
    parms[cn.PLT_YLABEL] = makeLabel(other_transfer, 
        self._mutation_column)
    parms[cn.PLT_XLABEL] = ""
    xpos = 0.8*len(mutations)
    ypos = 0.05*len(mutations)
    ax.text(xpos, ypos, makeLabel(transfer, 
        self._mutation_column))
    plot = ax.pcolor(df_plot, cmap='jet', vmin=COLORBAR_MIN,
        vmax=COLORBAR_MAX)
    labels = df_plot.columns.tolist()
    if parms.isTrue(cn.PLT_XAXISTICKTOP):
      ax.xaxis.tick_top()
    ax.set_xticks(np.arange(0.5, len(labels)))
    ax.set_xticklabels(labels, rotation=90)
    ax.set_yticks(np.arange(0.5, len(mutations)))
    ax.set_yticklabels(mutations)
    if parms.isTrue(cn.PLT_COLORBAR):
      fig.colorbar(plot, cmap='jet')
    parms.do(is_plot=is_plot)
    return df_plot
