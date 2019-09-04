"""Constructs data used to calculate the fractional co-occurrence of mutations."""

from microbepy.common import constants as cn
from microbepy.common.isolate import Isolate
from microbepy.common import util
from microbepy.correlation import genome_correlation
from microbepy.data.model_data_provider import ModelDataProvider
from microbepy.data import util_data

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


MIN_FRACTION = 0.25
THRESHOLD_FRAC = 0.2

# The following lines are excluded because they have
# only 1 transfer
EXCLUDED_LINES = ["HE2", "HR1", "UA2", "UE2"]


########################################################################
class MutationCofraction(object):

  def __init__(self, mutation_column=cn.GGENE_ID, species=None):
    """
    :param str mutation_column:
    """
    self._mutation_column = mutation_column
    self._species = species
    self.lines = self.getLines()
    self.ordered_mutations = self._getOrderedMutations()
    self.transfers = self._getTransfers()

  def _getTransfers(self):
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

  def getLines(self):
    """
    Obtains the lines present for one or both species.
    :return list-str: lines sorted by name
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
    if self._species is not None:
      result = query(species)
    else:
      result_dvh = query(cn.SPECIES_MIX_DVH)
      result_mmp = query(cn.SPECIES_MIX_MMP)
      result = list(set(result_dvh).intersection(result_mmp))
    result = [l for l in result if not l in EXCLUDED_LINES]
    result.sort()
    return result

  def makeLineDF(self, permitted_mutations=None,
       transfer=cn.TRANSFER_DEFAULT):
    """
    :param list-str permitted_mutations:
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
    # Eliminate excluded lines
    constraint = lambda r: not r[cn.LINE] in EXCLUDED_LINES
    df = util.selectRows(df, [constraint])
    # Construct the matrix with all relevant lines
    df_matrix = util.makeMatrix(df, row_name=self._mutation_column,
        column_name=cn.LINE, value_name=cn.FREQ)
    lines = self.getLines()
    lines = [l for l in lines if not l in EXCLUDED_LINES]
    for line in set(lines).difference(df_matrix.columns):
      df_matrix[line] = 0.0
    lines.sort()
    df_matrix = df_matrix[lines]
    # Select the species
    if self._species is not None:
      df_transpose = df_matrix.transpose()
      self.__class__._pruneColumns(df_transpose, self._species)
      df_matrix = df_transpose.transpose()
    # Adjust for permitted mutations
    if permitted_mutations is None:
      permitted_mutations = set(self._getFrequentMutations())
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
    df_matrix = df_matrix.sort_index()
    #
    return df_matrix

  def _getFrequentMutations(self, min_lines=2):
    """
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
    if self._species is not None:
      mutations = [m for m in mutations if m[0] == self._species]
    mutations.sort()
    return mutations

  def makeCofractionDF(self, transfer=cn.TRANSFER_DEFAULT,
      threshold_frac=THRESHOLD_FRAC,
      is_difference_frac=False,
      other_transfer=None):
    """
    Constructs a dataframe of the fraction of lines in which
    pairs of mutations occur.
    :param float theshold_frac: threshold fraction for co-occurence
    :param int transfer:
    :param bool is_difference_frac: Multiple other_transfer
        times 1 - its value in transfer.
    :return pd.DataFrame: columns and index are mutations
        values are fraction of lines in which mutations co-occur
    """
    def makeDF(transfer):
      df_line = self.makeLineDF(transfer=transfer)
      df_binary = df_line.applymap(
          lambda v: 0 if np.isnan(v) else v)
      df_binary = df_line.applymap(
          lambda v: 1.0 if v > threshold_frac else 0)
      if df_binary.isnull().sum().sum() > 0:
        raise ValueError("Unexpected null or nan.")
      return df_binary.transpose()
    #
    if other_transfer is None:
      other_transfer = transfer
    #
    df_binary_transfer = makeDF(transfer)
    df_binary_other = makeDF(other_transfer)
    col_transfer = df_binary_transfer.columns
    col_other = df_binary_other.columns
    df_binary_transfer = util_data.addRowsColumns(df_binary_transfer,
        col_other, 0, is_columns=True, is_rows=True)
    df_binary_other = util_data.addRowsColumns(df_binary_other,
        col_transfer, 0, is_columns=True, is_rows=True)
    if is_difference_frac:
      df_binary_other = df_binary_other * (1 - df_binary_transfer)
    # 
    df_counts = df_binary_transfer.T.dot(df_binary_other)
    df_result = df_counts.applymap(lambda v: v / len(self.lines))
    return df_result

  def _getOrderedMutations(self):
    """
    Calculate the mutation ordering from clustermap.
    This is done relative to 1K generations mutations.
    :param pd.DataFrame df_plot: mutations are rows and columns
    Notes:
      1. Closes current plot
    """
    df = self.makeCofractionDF(transfer=cn.TRANSFER_1000G,
        threshold_frac=0)
    df = df.fillna(0)
    cg = sns.clustermap(df)
    plt.close()
    return [df.index[i] for i in cg.dendrogram_row.reordered_ind]

  def makeCofractionDifferencedDF(self,
      transfer=cn.TRANSFER_DEFAULT,
      threshold_frac=THRESHOLD_FRAC,
      other_transfer=None):
    """
    Calculates the threshold difference of fractions between
    two transfer times.
    """
    def makeDF(tfr):
      df = self.makeCofractionDF(transfer=tfr,
          other_transfer=tfr)
      return util_data.addRowsColumns(df, self.ordered_mutations, 0,
          is_columns=True, is_rows=True)
    #
    df_this = makeDF(transfer)
    df_other = makeDF(other_transfer)
    df = df_this - df_other
    return df
