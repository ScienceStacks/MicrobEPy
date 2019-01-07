""" Transforms predictor variables.  """


import microbepy_init
import constants as cn
import util

import numpy as np
import pandas as pd

MIN_CORR = 1  # Minimum correlation for two mutations to be
                 # considered correlated and therefore grouped
MAX_FRAC = 0.90  # Maximum faction of cultures in which a mutation
                # can be present
MIN_DEG_FREE = 10  # Minimum degrees of freedom after parameter
                   # estimation and cross validation
MIN_NUM_LINES = 6
#


class PredictorTransformer(object):
  """
  Transforms a predictor variable matrix, df_X.
  The resulting predictor matrix is self.df_X.
  """

  df_culture_isolate = None

  def __init__(self, df_X, mutation_column):
    """
    :param pd.DataFrame df_X: indexed by cn.KEY_CULTURE,
        columns are predictor variables.
    :param str mutation_column: column containing the mutation
    """
    cls = self.__class__
    self.df_X = df_X.copy()
    self._mutation_column = mutation_column
    if cls.df_culture_isolate is None:
      cls._makeCultureIsolate()

  ####### Transformer Methods: Workflows #######

  def transformLowFrequencyColumns(self, max_columns):
    """
    :param int min_isolate_pairs:
    """
    self.makeMutationGroups()
    self.filterHighFrequencyMutations()
    self.filterLowFrequencyColumns(max_columns)

  def transformOnlyLowFrequencyColumns(self, max_columns):
    """
    :param int min_isolate_pairs:
    """
    self.filterLowFrequencyColumns(max_columns, min_shared_isolates=4)

  def transformDefault(self, **kwargs):
    """
    Groups mutations by their correlation and removes very frequent mutations.
    """
    self.filterForMultipleLines(**kwargs)
    self.makeMutationGroups()
    self.filterHighFrequencyMutations()

  ####### Filter Methods: Eliminate columns #######

  def filterForMultipleLines(self, min_num_lines=MIN_NUM_LINES, **kwargs):
    """
    Exclude mutations that are not in multiple lines.
    """
    mutations = self._getLowLineOccurrenceMutations(min_num_lines)
    for col in self.df_X:
      if col in mutations:
        del self.df_X[col]

  @classmethod
  def _makeCultureIsolate(cls):
    """
    Creates mapping between culture and isolate pairs
    :return pd.DataFrame:
      cn.KEY_CULTURE, cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP
    """
    query = '''
    select distinct key_culture,
            key_isolate as key_isolate_dvh, 
            key_isolate_mmp 
          from genotype_phenotype,
            (select distinct key_isolate as key_isolate_mmp, 
                line as line_mmp,
                key_culture as key_culture_mmp from genotype_phenotype 
              where species='M' 
                and line_mmp != 'AN') sub 
          where species_mix = 'B' 
              and species = 'D' 
              and key_culture_mmp = key_culture
              and line != 'AN'
    '''
    cls.df_culture_isolate = util.readSQL(query)

  def _getLowLineOccurrenceMutations(self, min_num_lines):
    """
    Finds the mutations that occur infrequently in lines.
    :param int min_num_lines:
    :return list-of-str:
    """
    query = '''
      select %s, cnt
        from (select distinct %s, count(distinct line) as cnt from genotype
         where is_an_mutation = 0
         group by %s) sub
        where cnt < %d
        order by cnt
      ''' % (self._mutation_column, self._mutation_column, 
          self._mutation_column, min_num_lines)
    df = util.readSQL(query)
    mutations = df[self._mutation_column].tolist()
    return mutations
 
  def makeMutationGroups(self, min_corr=MIN_CORR):
    """
    Creates groups of mutations, adjusting the predictor
    :param float min_corr: Minimum correlation to be considered for grouping
    """
    # Group together correlated columns in X
    def aggregateFunc(df):
      """
      Function that aggregates the columns of a dataframe 
      This aggregation does a boolean OR
      for binary valued columns.
      """
      result = df[df.columns[0]]  # Initialization
      for col in df.columns:
        result = result + df[col]
      result = result.apply(lambda v: min(v, 1))
      return result
    # Form groups recursively
    for _ in range(10):
      classes = util.findCorrelatedColumns(self.df_X, min_corr=min_corr)
      if len(classes) == 0:
        break
      for group in classes:
        util.aggregateColumns(self.df_X, group, aggregateFunc)

  def filterHighFrequencyMutations(self, max_frac=MAX_FRAC):
    """
    Eliminates mutation that occur too frequently since they
    have little regression leverage.
    :param float max_frac: Maximum fraction for the mutation
    """
    max_count = int(len(self.df_X)*max_frac)
    df_count = self.df_X.sum()
    columns = df_count[df_count >= max_count].index.tolist()
    for col in columns:
      del self.df_X[col]

  def filterLowFrequencyColumns(self, max_columns,
      min_shared_isolates=2):
    """
    Filters columns based on occurrences in isolate pairs.
    (aggregation of rows)
    :param int max_columns:
    """
    cls = self.__class__
    #
    df = self.df_X.copy()
    df = df.reset_index()
    if not set([cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP]).issubset(
        df.columns):
      df = df.merge(cls.df_culture_isolate, on=cn.KEY_CULTURE,
            how='inner')
    df_X = df.groupby([cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP]).mean()
    df_sum = df_X.sum()
    df_sum = df_sum[df_sum > min_shared_isolates]
    df_sum = df_sum.sort_values(ascending=False)
    upper = min(max_columns, len(df_sum))
    columns = df_sum.index[0:upper]
    self.df_X = self.df_X[columns].copy()
