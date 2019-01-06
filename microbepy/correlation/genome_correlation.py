"""
Provides correlations beween genomes for different types of genome instances.

Both genome and genome instance have different levels of granularity. A
genome can be at the granularity of a mutation key, mutation position,
or a gene. A genome instance can be an isolate or a line.

Key concepts:
  instance name - name of a column in which an instance occurs
  categorical name - name of a categorical variable whose values become columns

Example of plotting a correlation heatmap:
  genome_correlation = GenomeCorrelation()
  # Create a correlation matrix
  df = genome_correlation.makeCorrelationDF(row_name=KEY_ISOLATE, column_name=KEY_MUTATION)
  genome_correlation.plotHeatmap(df)
"""

import matplotlib.pyplot as plt

import __init__
import constants as cn
import pandas as pd
import numpy as np
import util
from isolate import Isolate
import correlation_statistic as cs

import copy

PLOT_WIDTH = 20
PLOT_HEIGHT = 10
# Columns
MUTE1 = "mute1"
MUTE2 = "mute2"
CORR = "corr"
COUNT_MUTATION = "count_mutation"
COUNT_ISOLATE = "count_isolate"

# Column Values
ALL = "all"

# Fonts
FONTSIZE_TITLE = 20 
FONTSIZE_AXIS = 14

#
TICKS_DEFAULT = [-1, 1.25, 0.25]


########################################
# Functions
########################################
def makeBinaryDF(df, instance_name, categorical_name, constraints):
  """
  Creates a binary dataframe for categorical variables.
  Rows are instances of occurrences. Columns are values of the
  categorical variable that occur with an instance.
  :param pd.DataFrame df: has column for instance and the categorical variable.
  :param str instance_name: Name of the column that defines an instance
                            to a assess correlation
  :param str categorical_name: Name of the column whose categorical values
                               are correlated
  :param list-of-function constraints: list of boolean function
                                       on self.df_base
  :return pd.DataFrame: instance_name, covar_name
  """
  df_sub = util.selectRows(df, constraints)
  return util.makeMatrix(df_sub, row_name=instance_name, 
      column_name=categorical_name)

def makeSiglvlDF(df_binary, df_other=None):
  """
  Creates a DF with entries that are the significance level of
  co-occurrences of an attribute (a categorical value)
  :param pd.DataFrame df_binary: A binary matrix as a dataframe that
      has the columns as the correlation variable and rows as samples
      and cell values are boolean.
  :param pd.DataFrame df_other: An optional second binary matrix
      structured as the first. If present, this matrix represents
      the columns in the result.
  :return pd.DataFrame: rows and columns are sample attributes,
      cell values are significance level of co-occurrences.
  """
  if df_other is None:
    df_other = df_binary
  if not set(df_binary.columns) == set(df_other.columns):
    raise ValueError("Dataframes must have the same columns")
  if not set(df_binary.index) == set(df_other.index):
    raise ValueError("Dataframes must have the same index")
  #
  num = len(df_binary)
  # Select common mutations
  mutations = list(set(df_binary.columns).intersection(df_other.columns))
  # Compute the count of co-occurrences for each attribute pair
  df_transpose = df_binary.copy()[mutations]
  df_transpose = df_transpose.transpose()
  df_corr_count = df_transpose.dot(df_other[mutations])
  # Compute the significance levels
  df_count_row = df_binary.sum(axis=0)
  df_count_column = df_other.sum(axis=0)
  items = df_binary.columns.tolist()
  rows = []
  for item1, cur_row in df_corr_count.iterrows():
    row = cur_row.copy(deep=True)
    row = row.astype(float)
    n1 = int(df_count_row[item1])
    for item2 in cur_row.index:
      n2 = int(df_count_column[item2])
      k = int(df_corr_count.loc[item1, item2])
      row[item2] = cs.calcCumlProb(num, n1, n2, k)
    rows.append(row)
  df_result = pd.DataFrame(rows)
  return df_result


########################################
# Classes
########################################
class GenomeCorrelation(object):

  """
  The base dataframe have the following columns:
     cn.COUNT - set to 1
     cn.KEY_MUTATION 
     cn.GENE_ID
     cn.GENE_POSITION
     cn.KEY_ISOLATE
     cn.COMMUNITY - community represented by the isolate
     cn.SPECIES
     cn.LINE
  """

  def __init__(self,
      plot_width=PLOT_WIDTH,
      plot_height=PLOT_HEIGHT,
      instance_name=cn.KEY_ISOLATE,
      categorical_name=cn.GGENE_ID,
      constraints=None,
      is_siglvl=True,
      is_test=False
      ):
    """
    :param str sequence: Data column for which sequence 
                         is considered
    :param int plot_width:
    :param int plot_height:
    :param str instance_name:
    :param str categorical_name:
    :param list-of-function contraints: boolean function on rows
        of the base dataframe
    :param bool is_siglvl: Use significance level in tests
        else use correlation
    :param bool is_test: True if this is a test invocation
    """
    self.plot_width = plot_width
    self.plot_height = plot_height
    self.instance_name = instance_name
    self.categorical_name = categorical_name
    self.is_test = is_test
    self._is_siglvl = is_siglvl
    self.df_base = self.__class__.makeBaseDF()
    self.df_binary = makeBinaryDF(self.df_base, self.instance_name, 
        self.categorical_name, constraints)
    self.df_binary.drop_duplicates(inplace=True)

  @classmethod
  def makeBaseDF(cls):
    """
    Creates a dataframe with the columns
      cn.KEY_MUTATION, cn.GENE, cn.GGENE_ID,
      cn.KEY_ISOLATE, cn.COUNT, cn.GENE_POSITION,
      cn.COMMUNITY, cn.SPECIES, cn.LINE
    Notes:
      Excludes ancestral mutations
    """
    query = '''
    select distinct key_mutation, gene_id, ggene_id, key_isolate, 
        species, line, gene_position
      from genotype 
      where species is not null
        and key_mutation is not null
        and is_an_mutation = 0
        and is_low_coverage_isolate = 0
    '''
    df_result = util.readSQL(query)
    df_result[cn.COUNT] = 1
    df_result[cn.COMMUNITY] = [Isolate.create(s).getCommunity()
        for s in df_result[cn.KEY_ISOLATE]]
    return df_result
      
  def makeCorrelationDF(self):
    """
    Creates a correlation dataframe.
    :return pd.DataFrame: columns are values of cn.KEY_MUTATION
    """
    if self._is_siglvl:
      df_result = makeSiglvlDF(self.df_binary)
    else:
      df_result = self.df_binary.corr()
      # Get nan values if a categorical_name instance is always present
      df_result = df_result.applymap(
          lambda v: 1.0 if np.isnan(v) else v)
    is_nan, columns = util.isNanInDataFrame(df_result)
    if is_nan:
      import pdb; pdb.set_trace()
    return df_result

  def plotMutationPairsHist(self, species=None):
    """ 
    Plots a histogram of the occurrence of mutation pairs in isolates.
    :param str species:
    :param bool self.is_test: doesn't do plot; returns constructed dataframe
    """
    if species is None:
      constraints = []
    else:
      constraints = [lambda r: r[MUTE1][0] == species]

    query = ''' 
    select distinct mute1, mute2, count(isolate1) as count from 
      (
        select distinct key_mutation as mute1, key_isolate as isolate1 from genotype
        where key_mutation is not null and key_isolate is not null and species is not null
        and is_an_mutation = 0 and transfer = 152 
        and is_low_coverage_isolate = 0
    ) sub1,
    (
        select distinct key_mutation as mute2, key_isolate as isolate2 from genotype
        where key_mutation is not null and key_isolate is not null and species is not null
        and is_an_mutation = 0 and transfer = 152
        and is_low_coverage_isolate = 0
    ) sub2
     where isolate1 = isolate2 and mute1 != mute2
     group by mute1, mute2
     order by count, mute1, mute2
     '''
    df = util.readSQL(query)
    df = util.selectRows(df, constraints)
    plt.hist(df[cn.COUNT], 30, density=True, cumulative=True, facecolor='g', alpha=0.75)
    plt.xlabel('# Isolates')
    plt.ylabel('Probability')
    if species is None:
      species = "%s and %s" % (cn.SPECIES_MIX_DVH, cn.SPECIES_MIX_MMP)
    title = "Distribution of number of isolates in which mutations co-occur for %s" % species 
    plt.title(title)
    upper = max(df[cn.COUNT].tolist())
    plt.axis([0, upper, 0, 1.0])
    plt.grid(True)
    if self.is_test:
      return df
    else:
      plt.show()

  def makeCountDF(self):
    """
    Creates a dataframe that counts mutations and isolates by line.
    :return pd.DataFrame: cn.LINE, COUNT_MUTATION, COUNT_ISOLATE
    """
    def isMatch(a, b):
      if (a == ALL) or (b == ALL):
        return True
      else:
        return a == b

    lines = self.df_base[cn.LINE].unique().tolist()
    speciess = self.df_base[cn.SPECIES].unique().tolist()
    speciess.append(ALL)
    lines.append(ALL)
    counts = {k: [] for k in [cn.LINE, cn.SPECIES, COUNT_MUTATION,
              COUNT_ISOLATE]}
    for line in lines:
      for species in speciess:
        count_keys = {}
        for key in [cn.KEY_MUTATION, cn.KEY_ISOLATE]:
          df = self.df_base[[cn.LINE, cn.SPECIES, key]]
          rows = [r for _,r in df.iterrows() 
              if (isMatch(r[cn.LINE], line)) and isMatch(r[cn.SPECIES], species)]
          count_keys[key] = len(set(pd.DataFrame(rows)[key]))
        counts[cn.LINE].append(line)
        counts[cn.SPECIES].append(species)
        counts[COUNT_MUTATION].append(count_keys[cn.KEY_MUTATION])
        counts[COUNT_ISOLATE].append(count_keys[cn.KEY_ISOLATE])
    df_result = pd.DataFrame(counts)
    df_result = df_result[[cn.LINE, cn.SPECIES,
        COUNT_MUTATION, COUNT_ISOLATE]]
    return df_result

  def plotGroupbyCountHist(self, groupby_column, count_column, species=None):
    """ 
    Plots a histogram of the count of occurrences.
    :param str groupby_column: either KEY_MUTATION or KEY_ISOLATE
    :param str column_column: either KEY_MUTATION or KEY_ISOLATE
    :param str species: only plot the species
    These are non-ancestral mutations.
    """
    ALL = 'all'
    SPECIES_NAME = {cn.SPECIES_MIX_DVH: 'DVH', cn.SPECIES_MIX_MMP: 'MMP'}
    COLUMN_NAME = {cn.KEY_MUTATION: "Mutations", cn.KEY_ISOLATE: "Isolates"}
    #
    groupby_name = str(COLUMN_NAME[groupby_column])
    groupby_name = groupby_name[:-1]  # Drop the 's'
#
    fig = plt.figure()
    fig.set_size_inches(PLOT_WIDTH, PLOT_HEIGHT)
    constraints = []
    for col in [groupby_column, count_column]:
      constraints.append(lambda r: not util.isNull(col))
    if species is not None:
      constraints.append(lambda r: r[cn.SPECIES] == species)
    lines = self.df_base[cn.LINE].unique().tolist()
    lines.insert(0, ALL)
    idx = 1
    plt.tight_layout()
    for line in lines:
      all_constraints = list(constraints)
      if line != ALL:
        all_constraints.append(lambda r: r[cn.LINE] == line)
      plt.subplot(2, 2, idx)
      df = util.selectRows(self.df_base, all_constraints)
      df = df[[groupby_column, count_column]]
      df.drop_duplicates(inplace=True)
      df_count = df.groupby(groupby_column).count()
      df_count.reset_index(inplace=True)
      plt.hist(df_count[count_column], 30, 
          density=True, cumulative=True, facecolor='g', alpha=0.75)
      if species is None:
        species = "%s and %s" % (cn.SPECIES_MIX_DVH, cn.SPECIES_MIX_MMP)
      if idx == 1:
        title = "%s per %s for %s by %s" % (
            COLUMN_NAME[count_column], 
            groupby_name,
            SPECIES_NAME[species], 
            line,
            )
      else:
        #num_groupby = len(df_count[groupby_column].unique())
        #title = "%s (%d)" % (line, num_groupby)
        title = "%s" % line
      if idx == 3:
        plt.ylabel("Culmulative %ss" % groupby_name)
        plt.xlabel('# %s' % COLUMN_NAME[count_column])
      else:
        plt.xlabel('')
        plt.ylabel('')
      plt.title(title)
      upper = max(df_count[count_column])
      plt.axis([0, upper, 0, 1.0])
      plt.grid(True)
      idx += 1
    if self.is_test:
      return df_count
    else:
      plt.show()

  def plotHeatmap(self, df, title="", ticks=None):
    """
    Constructs a heatmap for a matrix such as the output of
    makeCorrelationDF.
    :param pd.DataFrame df: Structured as a matrix for which
        a value is plotted for each combination of row and column
    :param str title:
    :param list-of-float ticks: how ticks are spaced - start, end, increment
    """
    if ticks is None:
      ticks = TICKS_DEFAULT
    if self._is_siglvl:
      # Compute one minus significance level
      df = df.copy()
      df = df.applymap(lambda x: 1 - x)
    # Create the figure
    fig, ax = plt.subplots()
    heatmap = ax.pcolor(df, cmap='seismic')
    fig.set_size_inches(self.plot_width, self.plot_height)
    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(df.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(df.shape[0]) + 0.5, minor=False)
    ax.set_title(title, fontsize=FONTSIZE_TITLE)
    # Define row and column labels
    column_labels = df.columns.tolist()
    row_labels = column_labels
    _ = ax.set_xticklabels(column_labels, minor=False, rotation=90)
    _ = ax.set_yticklabels(row_labels, minor=False)
    # Add the colorbar as a legend
    cbar = plt.colorbar(heatmap)
    cbar.set_ticks(np.arange(ticks[0], ticks[1], ticks[2]))
    if not self.is_test:
      plt.show()
