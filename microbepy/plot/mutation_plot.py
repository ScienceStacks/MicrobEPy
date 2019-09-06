"""Provides plots of  mutations for Isolates and Lines."""

from microbepy.common import constants as cn
from microbepy.common.dataframe_sorter import DataframeSorter
from microbepy.common.isolate import Isolate
from microbepy.common import util
from microbepy.correlation import genome_correlation
from microbepy.data.model_data_provider import ModelDataProvider
from microbepy.data import util_data
from microbepy.plot.mutation_cofraction import MutationCofraction
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
MAX_LINES = 9
MIN_FRACTION = 0.25
THRESHOLD_FRAC = 0.2
MAX_SIGLVL = 0.01
COLORBAR_MIN = 1.0
COLORBAR_MAX = 4.0


class MutationLinePlot(object):
  """
  Plot mutations by occurrences within Lines.
  """

  def __init__(self, mutation_column=cn.GGENE_ID, species=None,
      is_plot=True):
    """
    :param str mutation_column:
    :param bool is_plot:
    """
    self._mutation_column = mutation_column
    self._is_plot = is_plot
    self._species = species
    self.cofraction = MutationCofraction(species=self._species,
        mutation_column=mutation_column)

  def plotTransfers(self,
        parms=PlotParms(is_initialize=False),
        is_unit_fraction = False,
        is_cluster_mutations=True):
    """
    Does a stacked bar plot of mutation frequency for all transfers.
    :param bool is_unit_fraction: round fraction to 1
    :param bool is_cluster_mutations: Group similar mutations together
    :return pd.DataFrame: row=mutation, col=line + transfer, value is fraction
    """
    permitted_mutations = self.cofraction.ordered_mutations
    transfers = self.cofraction.transfers
    num_transfers = len(transfers)
    fig, axes = plt.subplots(nrows=num_transfers, ncols=1)
    dfs = []
    for idx, transfer in enumerate(transfers):
      parms[cn.PLT_YTICKLABELS] = True
      if self._species is None:
        parms[cn.PLT_TITLE] = "%d" % transfer
      else:
        parms[cn.PLT_TITLE] = "%s, %d" % (self._species, transfer)
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
      df = self.plotLine(transfer, 
          parms=parms, is_plot=False,
          ax=axes[idx], permitted_mutations=permitted_mutations,
          is_unit_fraction=is_unit_fraction)
      df[cn.TRANSFER] = transfer
      dfs.append(df)
    if self._is_plot:
      plt.show()
    return pd.concat(dfs)
    

  def plotLine(self, transfer, 
      parms=PlotParms(is_initialize=False),
      is_unit_fraction=False,
      is_plot=None, ax=None, permitted_mutations=None):
    """
    Does a stacked bar plot of mutation frequency by line
    with colors
    :params int transfer:
    :params PlotParms parms:
    :params Axis ax: axis to use in plot
    :param list-str permitted_mutations: to use and how they
       are ordered if None, then use alphabetical order
    :param bool is_unit_fraction: round non-zero fraction to 1
    :return pd.DataFrame: row=mutation, col=line, value is fraction
    """
    if is_plot is None:
      is_plot = self._is_plot
    parms.setTrueIfAbsent(cn.PLT_XLABEL)
    parms.setTrueIfAbsent(cn.PLT_XTICKLABELS)
    #
    df_plot = self.cofraction.makeLineDF(
        permitted_mutations=permitted_mutations,
        transfer=transfer)
    if is_unit_fraction:
      df_plot = df_plot.applymap(
          lambda v: 1 if v> MIN_FRACTION else v)
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
      title = "%s Mutations" % SPECIES[self._species]
    xpos = int(len(df_plot)*0.5)
    ypos = MAX_LINES - 3
    ax.text(xpos, ypos, title, fontsize=FONTSIZE_TITLE)
    ax.set_ylim([0, MAX_LINES])
    if parms.isTrue(cn.PLT_YLABEL):
      if is_unit_fraction:
        label = "No. Lines"
      else:
        label = "Fraction"
      ax.set_ylabel(label , fontsize=FONTSIZE_LABEL)
    if parms.isTrue(cn.PLT_XLABEL):
      ax.set_xlabel(self._mutation_column, fontsize=FONTSIZE_LABEL)
    if parms.isTrue(cn.PLT_LEGEND):
      ax.legend(loc=(1,2))
      #ax.legend()
    if is_plot:
      plt.show()
    return df_plot

  def _makeMutationSiglvlMatrix(self,
       transfer=cn.TRANSFER_DEFAULT, 
       other_transfer=None, min_fraction=MIN_FRACTION):
    """
    Creates a significance level matrix for mutations.
    :param int transfer: transfer time for row mutations
    :param int other_transfer: transfer time for column mutations
    :param float min_fraction: minimum fractional occurrence of
        a mutation within a line for it to be considered
    :return pd.DataFrame: row index and columns are mutations
    """
    def makeDF(transfer):
      df_line = self.cofraction.makeLineDF(transfer=transfer)
      df_binary = df_line.applymap(
          lambda v: 0 if np.isnan(v) else v)
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
      max_siglvl=MAX_SIGLVL):
    """
    Constructs a the dataframe used for heatmap.
    :param int transfer:
    :param float max_siglvl:
    :return pd.DataFrame: mutations, mutations,
        values are -log10 significance level
    """
    df_matrix = self._makeMutationSiglvlMatrix(transfer=transfer,
        other_transfer=other_transfer)
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

  def plotCofractions(self, is_time_lag=False,
      threshold_frac=THRESHOLD_FRAC,
      is_difference_frac=False,
      is_differenced=False,
      is_compress=False,
      parms=PlotParms(), **kwargs):
    """
    Does a subplots of the fraction of lines in which mutations co-occur.
    :param bool is_time_lag: construct time lag subplots
    :param bool is_differenced: Computes the difference in
        count fractions
    :param dict kwargs: non-transfer parameters passed to next level
    :return dict: key is pair of transfers, value is data_frame
    """
    def funcDF(transfer, other_transfer):
      if is_differenced:
        df = self.cofraction.makeCofractionDifferencedDF(
            transfer=transfer, other_transfer=other_transfer,
            threshold_frac=threshold_frac)
      else:
        df = self.cofraction.makeCofractionDF(transfer=transfer,
            is_difference_frac=is_difference_frac,
            other_transfer=other_transfer)
      if is_compress:
        df.dropna(axis=0, how='all', inplace=True)
        df.dropna(axis=1, how='all', inplace=True)
      return df
    #
    return self._plotTransfers(funcDF, is_time_lag, 
        parms=parms, heat_range=[0, 1.0], **kwargs)

  def plotSiglvls(self, is_time_lag=False, max_siglvl=MAX_SIGLVL,
      parms=PlotParms(), **kwargs):
    """
    Does a subplots of mutation correlation significance levels.
    :param bool is_time_lag: construct time lag subplots
    :param dict kwargs: non-transfer parameters passed to next level
    :return dict: key is pair of transfers, value is data_frame
    """
    def funcDF(transfer, other_transfer):
      return self._plotSiglvlDF(transfer=transfer,
          max_siglvl=max_siglvl,
          other_transfer=other_transfer)
    #
    return self._plotTransfers(funcDF, is_time_lag, 
        parms=parms, 
        heat_range = [COLORBAR_MIN, COLORBAR_MAX],
        **kwargs)

  def _plotTransfers(self, funcDF, is_time_lag, 
      parms=PlotParms(), **kwargs):
    """
    Does a subplots of mutation mutations over transfers.
    :param Function funcDF: has kwargs transfer, other_transfer;
        returns a dataframe of mutations as columns and index;
        values are used in the heatmap.
    :param bool is_time_lag: construct time lag subplots
    :param dict kwargs: non-transfer parameters passed to next level
    :return dict: key is pair of transfers, value is data_frame
    """
    NCOLS = 3
    plot_pos = {1:1, 2:3, 3:4, 4:6}
    NPLOTS = 6
    transfers = self.cofraction.transfers
    if is_time_lag:
      pairs = [p for p in zip(transfers[:-1], transfers[1:])]
    else:
      pairs = [p for p in zip(transfers[:-1], transfers[:-1])]
    #
    # Calculate the column order
    df = funcDF(transfer=cn.TRANSFER_1000G,
        other_transfer=cn.TRANSFER_1000G)
    df = df.fillna(0)
    # Set up for plots
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
      transfer = pair[0]
      other_transfer = pair[1]
      df = funcDF(transfer=transfer, other_transfer=other_transfer)
      df = df.applymap(lambda v: np.nan if v == 0 else v)
      self._plotTransferCompare(df, 
          transfer=transfer, other_transfer=other_transfer,
          ordered_columns=self.cofraction.ordered_mutations,
          is_center_colorbar=True,
          fig=fig, ax=ax, parms=parms, is_plot=is_plot, **kwargs)
      result[pair] = df
    return result

  def plotSiglvl(self, max_siglvl=MAX_SIGLVL, 
      transfer=cn.TRANSFER_DEFAULT,
      other_transfer=None,
      is_center_colorbar = True,
      **kwargs):
    """
    Constructs a heatmap of the mutation coocurrence significance
    levels.
    :param float max_siglvl: maximum significance level
    :return pd.DataFrame: columns, rows are mutations
    """
    df_plot = self._plotSiglvlDF(transfer=transfer,
        other_transfer=other_transfer,
        max_siglvl=max_siglvl)
    self._plotTransferCompare(df_plot, 
        heat_range = [COLORBAR_MIN, COLORBAR_MAX],
        ordered_mutations=self.cofraction.ordered_mutations,
        transfer=transfer, other_transfer=other_transfer,
        is_center_colorbar=is_center_colorbar,
        **kwargs)
    return df_plot

  def plotCofraction(self,
      threshold_frac=THRESHOLD_FRAC,
      transfer=cn.TRANSFER_DEFAULT,
      other_transfer=None,
      is_difference_frac=False,
      is_differenced=False,
      is_center_colorbar=True,
      is_compress=False,
      parms=PlotParms(),
      **kwargs):
    """
    Constructs a heatmap of the mutation coocurrence fractions.
    :param int transfer: Transfer for which plot is done
    :param bool is_differenced: Computes the difference in
        count fractions
    :param bool is_compress: Eliminate rows/columns
       with 0 values
    :return pd.DataFrame: columns, rows are mutations
    """
    if is_differenced:
      df = self.cofraction.makeCofractionDifferencedDF(
          threshold_frac=threshold_frac,
          transfer=transfer, other_transfer=other_transfer,
          **kwargs)
      df = df.applymap(lambda v: np.nan 
          if np.abs(v) < threshold_frac else v)
    else:
      df = self.cofraction.makeCofractionDF(transfer=transfer,
          is_difference_frac=is_difference_frac,
          other_transfer=other_transfer, **kwargs)
      df = df.applymap(lambda v: np.nan if v < threshold_frac else v)
    if is_compress:
      df.dropna(axis=0, how='all', inplace=True)
      df.dropna(axis=1, how='all', inplace=True)
      is_include_missing_mutations = False
    else:
      is_include_missing_mutations = True
    ordered_columns = self.cofraction.ordered_mutations
    self._plotTransferCompare(df, 
        heat_range=[0, 1.0],
        ordered_columns=ordered_columns,
        parms=parms,
        transfer=transfer, other_transfer=other_transfer,
        is_center_colorbar=is_center_colorbar,
        is_include_missing_mutations=is_include_missing_mutations,
        **kwargs)
    return df

  def _plotTransferCompare(self, 
      df_plot,
      heat_range,
      ordered_columns=None,
      is_center_colorbar=True,
      transfer=cn.TRANSFER_DEFAULT,
      other_transfer=None,
      ax=None,
      fig=None,
      is_include_missing_mutations=True,
      parms=PlotParms(),
      is_plot=None):
    """
    Constructs a heatmap comparing values for mutations from two transfers.
    :param pd.DataFrame df_plot: index and columns are mutations;
        values are plotted on the heatmap
    :param list-str ordered_columns: order in which columns appear
    :param bool is_center_colorbar: center the colorbar in the plot
    :param float, float: values on the heatmap range
    :param int transfer:
    :param int other_transfer: Allow comparisons across time
    :param Matplotlib.Axes ax:
    :param PlotParms parms: Parameters for the plot
    :param bool is_plot: Overrides constructor plotting directive
    :param bool is_include_missing_mutations:
    """
    def makeLabel(transfer, column, is_include_column=False):
      if is_include_column:
        label = "%d-%s" % (transfer, column)
      else:
        label = "%d" % transfer
      return label
    def setValue(a_dict, key, default):
      if not key in a_dict.keys():
        a_dict[key] = default
    #
    if is_plot is None:
      is_plot = self._is_plot
    elif not self._is_plot:
      is_plot = self._is_plot
    if ordered_columns is None:
      ordered_columns = list(set(df_plot.columns.tolist()).union(
          df_plot.index))
    # Do the plot
    if not cn.PLT_COLORBAR in parms:
      parms[cn.PLT_COLORBAR] = True
    if other_transfer is None:
        other_transfer = transfer
    if ax is None:
      if fig is None:
        fig = plt.figure(figsize=parms[cn.PLT_FIGSIZE])
      ax = fig.add_subplot(1, 1, 1)
    # Order the columns
    if is_include_missing_mutations:
      columns = df_plot.columns.tolist()
      missing_columns = set(ordered_columns).difference(columns)
      extended_ordered_columns = list(ordered_columns)
      extended_ordered_columns.extend(
          set(columns).difference(ordered_columns))
      for col in missing_columns:
        df_plot[col] = np.nan
        df_plot.loc[col, :] = np.nan
      df_plot = df_plot.reindex(extended_ordered_columns)
      df_plot = df_plot[extended_ordered_columns]
      rows = df_plot.columns.tolist()
      columns = df_plot.columns.tolist()
    else:
      extended_ordered_columns = ordered_columns
      rows = df_plot.index.tolist()
      columns = df_plot.columns.tolist()
    mutations = df_plot.columns.tolist()
    # Set up plot information
    parms[cn.PLT_XLABEL] = ""
    setValue(parms, cn.PLT_COLORBAR, True)
    xpos = 1.05*len(columns)
    ypos = -0.05*len(rows)
    parms[cn.PLT_XLABEL] = ""
    xlabel = makeLabel(other_transfer, self._mutation_column)
    parms[cn.PLT_YLABEL] = makeLabel(
        transfer, self._mutation_column)
    ax.text(xpos, ypos, xlabel, fontsize=parms.fontsize_label)
    #
    # Construct the plot
    plot = ax.pcolor(df_plot, cmap='jet', vmin=heat_range[0],
        vmax=heat_range[1])
    if parms.isTrue(cn.PLT_COLORBAR):
      if is_center_colorbar:
        # Colorbar positions: left, bottom, width, height
        cbaxes = fig.add_axes([.45, 0.2, 0.01, 0.5]) 
        cb = fig.colorbar(plot, cax = cbaxes, cmap='jet')
        cb.ax.tick_params(labelsize=parms.fontsize_label)
      else:
        cb = fig.colorbar(plot, cmap='jet')
        cb.ax.tick_params(labelsize=parms.fontsize_label)
    row_labels = df_plot.columns.tolist()
    col_labels = df_plot.index.tolist()
    if parms.isTrue(cn.PLT_XAXISTICKTOP):
      ax.xaxis.tick_top()
    ax.set_xticks(np.arange(0.5, len(row_labels)))
    ax.set_xticklabels(row_labels, rotation=90,
        fontsize=parms.fontsize_label)
    ax.set_yticks(np.arange(0.5, len(col_labels)))
    ax.set_yticklabels(col_labels,
        fontsize=parms.fontsize_label)
    #parms[cn.PLT_YLABEL] = ""
    parms.do(is_plot=False)
    if is_plot:
      parms[cn.PLT_YLABEL] = ""
      parms.do(is_plot=False)
      ylabel = makeLabel(transfer, self._mutation_column)
      xpos = -3
      ypos = 0.5*len(rows)
      ypos = -1
      ax.set_ylabel(ylabel, fontsize=parms.fontsize_label,
          x=xpos, y=ypos)
      #plt.show()
      parms.do(is_plot=is_plot)
    else:
      parms.do(is_plot=is_plot)
