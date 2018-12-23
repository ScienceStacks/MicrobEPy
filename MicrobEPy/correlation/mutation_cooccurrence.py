"""
Analyzes the co-occurrence of mutations in isolates.

Produces Combination Statistics - CoStatistic. 
A combination is a group
of Isolates. Combination statistics are statistics about
a combination that summarize
  - mutations common to the combination
  - the average value of dependent variable(s)
  - the range of values of dependent variable(s)

The CoStatistic dataframes are organized by keys, each key
is associated with a different dataframe:
    cn.MIN
       cn.ISOLATES - isolates in group
       cn.MUTATIONS - list of mutations in common
       cn.COUNT_MUTATIONS - number of mutations in common
       cn.COUNT_ISOLATES - number of isolates
       cn.MIN - minimum of the rate and yield
       cn.MAX - maximum of the rate and yield
    cn.RATE, cn.YIELD
       cn.ISOLATES - isolates in group
       cn.MUTATIONS - list of mutations in common
       cn.COUNT_MUTATIONS - number of mutations in common
       cn.COUNT_ISOLATES - number of isolates
       cn.AVG - average value in units of standard deviation
       cn.RNG - range in units of standard deviation

The CoStatistic dataframes are constructed from Sample dataframes
that are indexed by isolate pair and have the following columns:
    cn.GROUP, cn.RATE, cn.YIELD, mutations (with binary values)
"""

import constants as cn
from model_data_provider import ModelDataDualProvider,  \
    ModelDataProvider
from study_context import nextStudyContext
from range_constraint import RangeConstraint
import util

from collections import namedtuple
import itertools
from sklearn import linear_model
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

NUM_REPLICATIONS = int(1e4)
XAXIS_VALUES = [cn.RATE, cn.YIELD, cn.MIN, cn.MAX]



class CoStatistic(object):
  # Container for Cooccurrence Statistics
  SCHEMA = {
      cn.MIN: [cn.ISOLATES, cn.MUTATIONS, 
      cn.COUNT_MUTATIONS, cn.MIN, cn.MAX],
      cn.RATE: [cn.ISOLATES, cn.MUTATIONS, 
      cn.COUNT_MUTATIONS, cn.AVG, cn.RNG],
      cn.YIELD: [cn.ISOLATES, cn.MUTATIONS, 
      cn.COUNT_MUTATIONS, cn.AVG, cn.RNG],
      }
  ATTRIBUTES = SCHEMA.keys()
 
  def __init__(self):
    self.dfs = {}
    for key in self.__class__.SCHEMA.keys():
      self.dfs[key] = pd.DataFrame()

  def get(self, name):
    return self.dfs[name]

  def set(self, name, df):
    self.dfs[name] = df

  def values(self):
    return self.dfs.values()

  def concat(self, other):
    """
    The current CoStatistic is concatenated with another CoStatistic.
    :param CoStatistic other:
    """
    for key in self.__class__.SCHEMA.keys():
      self.dfs[key] = pd.concat([self.dfs[key], other.dfs[key]])


RegressionResult = namedtuple('RegressionResult', 
    ['predictions', 'rsq', 'slope'])

class MutationCooccurrence(object):

  def __init__(self, mutation_column=cn.GGENE_ID,
      provider=None, constraints=None, is_plot=True):
    """
    :param ModelDataDualProvider provider: if specified, has invoked do()
    :param bool is_plot: plots if True
    """
    self._mutation_column = mutation_column
    self._is_plot = is_plot
    if provider is None:
      provider = ModelDataDualProvider(self._mutation_column, 
          constraints=constraints)
      provider.do(
          transform_type=cn.TRANSFORM_ONLY_LOW_FREQUENCY_ISOLATES)
    self.df_X = provider.df_X
    self.df_ys = provider.df_ys
    self.df_y_stds = provider.df_y_stds
    self.isolate_dict = ModelDataProvider.getIsolatesFromIndices(
        self.df_X.index)

  @staticmethod
  def _combineIsolates(isolate_dict):
    isolates = []
    [isolates.extend(v) for v in isolate_dict.values()]
    return isolates

  def findWithRangeConstraints(self, rc_vector):
    """
    Finds the set of mutations common to isolates satisfying
    constraints on rate and yield.
    :param RangeConstraintVector rc_vector:
        has keys cn.RATE, cn.YIELD
    :return dict, list-str: isolate_dict, mutations
        dictionary keyed by cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP
    """
    cls = self.__class__
    #
    # Select rows that abide by the RangeConstraints
    df = pd.DataFrame({
        cn.RATE: self.df_ys[cn.RATE][cn.VALUE],
        cn.YIELD: self.df_ys[cn.YIELD][cn.VALUE],
        })
    indices = rc_vector.findSatisfiedRows(df)
    isolate_dict = ModelDataProvider.getIsolatesFromIndices(indices)
    superset_isolates = cls._combineIsolates(isolate_dict)
    #
    isolates, mutations = self.find(
        superset_isolates=superset_isolates)
    return isolates, mutations
        
  def find(self,
      superset_isolates=None, superset_mutations=None):
    """
    Finds co-occurring mutations for a collection of isolates
    :param list-str superset_isolates: set of isolates considered
    :param list-str superset_mutations: set of mutations considered
    :return list-str, list-str: 
        isolates satisfying the constraints that are present in data
        list of mutations shared by isolates
    Notes:
      1. Includes the mutations of paired isolates
    """
    cls = self.__class__
    #
    if superset_isolates is None:
      superset_isolates = self.isolate_dict[cn.KEY_ISOLATE_DVH]
      superset_isolates.extend(self.isolate_dict[cn.KEY_ISOLATE_MMP])
    if superset_mutations is None:
      superset_mutations = self.df_X.columns.tolist()
    # Select the rows
    sel = [(i[0] in superset_isolates) or (i[1] in superset_isolates)
        for i in self.df_X.index]
    df_X = self.df_X[sel]
    isolate_dict = ModelDataProvider.getIsolatesFromIndices(
        df_X.index)
    isolates = cls._combineIsolates(isolate_dict)
    # Select the columns
    columns = [c for c in df_X.columns if c in superset_mutations]
    df_X_final = df_X[columns]
    # FInd the common mutations
    ser = df_X_final.product()
    mutations = [i for i, r in ser.items() if ser[i] == 1]
    #
    return isolates, mutations

  def _makeCoStatisticFromSampleDF(self, df_sample):
    """
    Computes statistics for aggregations of isolate groups.
    :param pd.DataFrame df_sample: indexed by isolates 
       cn.GROUP, Mutation Columns, cn.RATE, cn.YIELD
    :return CoStatistic:
    """
    #
    df_group = df_sample.groupby(cn.GROUP)
    groups = df_group.groups
    #
    df_cooccur = df_group.prod()
    for col in [cn.RATE, cn.YIELD]:
        del df_cooccur[col]
    ser_cooccur = df_cooccur.sum(axis=1)  # Indicates presence of mutation
    # Find the common mutations for each group
    mutation_stgs = []
    for group in groups.keys():
      df = df_cooccur[df_cooccur.index == group]
      mutations = [m for m in df_cooccur.columns 
                  if df.loc[group, m] == 1]
      mutations.sort()
      mutation_stgs.append(str(mutations))
    # Calculate the isolate strings
    isolates = []
    for key, values in groups.items():
      size = len(groups[key])
      isolate_stg =  [str(v) for v in values]
      isolate_stg.sort
      isolates.append(str(isolate_stg))
    # Compute common data
    df_min = df_sample.groupby(cn.GROUP).min()
    df_max = df_sample.groupby(cn.GROUP).max()
    df_avg = df_sample.groupby(cn.GROUP).mean()
    # Compute RATE, YIELD
    def makeDF(depvar):
      return pd.DataFrame({
          cn.COUNT_MUTATIONS: ser_cooccur,
          cn.MUTATIONS: mutation_stgs,
          cn.ISOLATES: isolates,
          cn.AVG: df_avg[depvar],
          cn.RNG: df_max[depvar] - df_min[depvar],
          })
    #
    df_result_rate = makeDF(cn.RATE)
    df_result_yield = makeDF(cn.YIELD)
    # Compute cn.MIN dataframe
    df_result_min = pd.DataFrame({
        cn.COUNT_MUTATIONS: ser_cooccur,
        cn.MUTATIONS: mutation_stgs,
        cn.ISOLATES: isolates,
        })
    df_result_min[cn.MAX] = pd.concat(
        [df_result_rate[cn.RNG], df_result_yield[cn.RNG]],
        axis=1).max(axis=1)
    df_result_min[cn.MIN] = pd.concat(
        [df_result_rate[cn.RNG], df_result_yield[cn.RNG]],
        axis=1).min(axis=1)
    result = CoStatistic()
    result.set(cn.RATE, df_result_rate)
    result.set(cn.YIELD, df_result_yield)
    result.set(cn.MIN, df_result_min)
    # Add the count of isolates
    for df in result.values():
      df[cn.COUNT_ISOLATES] = df[cn.ISOLATES].apply(
          lambda v: len(eval(v)))
    #
    return result  

  def makeCoStatistic(self, size, is_resample=False, **kwargs):
    """
    Chooses sets of isolates of the specified size.
    Provides the distribution of the number of mutations in
    common in the randomly chosen sets.
    :param int size: set size
    :param dict **kwargs: arguments for resample
    :param bool is_resample: always use resampling
    :return CoStatistic:
    """
    if is_resample or (len(self.df_X)**size > 1e6):
      df_sample = self._makeResampleSampleDF(size, 
          **kwargs)
    else:
      df_sample = self._makeExplicitSampleDF(size)
    return self._makeCoStatisticFromSampleDF(df_sample)

  def _makeExplicitSampleDF(self, size):
    """
    Constructs a dataframe of samples by explicitly constructing
    the possible sets of a desired size.
    :param int size: set size
    :return pd.DataFrame: Sample DataFrame
    """
    df_X = self.df_X.copy()
    isolates = df_X.index.tolist()
    isolate_indices = range(len(isolates))
    df_X[cn.INDEX] = isolate_indices
    df_X[cn.YIELD] = self.df_ys[cn.YIELD][cn.VALUE]
    df_X[cn.RATE] = self.df_ys[cn.RATE][cn.VALUE]
    df_X = df_X.reset_index()
    df_X[cn.INDEX] = isolate_indices
    combination_iterator = itertools.combinations(isolate_indices, 
        size)
    indices = []
    groups = []
    group_num = 0
    for combination in combination_iterator:
      indices.extend(list(combination))
      groups.extend([group_num] * size)
      group_num += 1
    df_index = pd.DataFrame({
        cn.GROUP: groups,
        cn.INDEX: indices,
        })
    df_result = df_index.merge(df_X, on=cn.INDEX, how='inner')
    df_result = df_result.sort_values(cn.GROUP)
    tuples = [(r[cn.KEY_ISOLATE_DVH], r[cn.KEY_ISOLATE_MMP]) 
        for _,r in df_result.iterrows()]
    del df_result[cn.INDEX]
    df_result.index = tuples
    return df_result
    
  def _makeResampleSampleDF(self, size, 
      num_replications=NUM_REPLICATIONS):
    """
    Constructs a dataframe of samples by resampling.
    :param int size: set size
    :param int num_replications: number of replications
    :return pd.DataFrame: Sample DataFrame
    """
    # Construct the replicated dataframe
    df_base = self.df_X.copy()
    length = len(df_base)
    df_base[cn.RATE] = self.df_ys[cn.RATE][cn.VALUE]
    df_base[cn.YIELD] = self.df_ys[cn.YIELD][cn.VALUE]
    df_sample = pd.concat([df_base] * num_replications)
    groups = []
    [groups.extend([n] * length) for n in range(num_replications)]
    df_sample[cn.GROUP] = groups
    # Add the sort number
    df_sample[cn.SORT] = np.random.uniform(
        0, 1, length*num_replications) + df_sample[cn.GROUP]
    df_sample = df_sample.sort_values(cn.SORT)
    del df_sample[cn.SORT]
    # Take the first elements of each group
    sel_base = [[True] * size, [False] * (length - size)]
    sel_lists = sel_base * num_replications
    sel = []
    [sel.extend(v) for v in sel_lists]
    return df_sample[sel]

  def plot(self, size,
      columns=[cn.RATE, cn.YIELD, cn.MIN, cn.MAX],
      title=""):
    """
    Scatter plot that is overlayed with a regression line.
    :param int size: Set size to plot
    :param list-str columns: Columns in co-occurrence dataframe
        that are to be plotted.
    :param str title: plot title
    :return pd.DataFrame:
        cn.RSQ, cn.SLOPE, cn.VALUE (size), cn.XAXIS
    """
    def regress(x, y):
      """
      :return RegressionResult:
      """
      X = pd.DataFrame({
          'ones': [1] * len(x),
          'x': x
          })
      lr = linear_model.LinearRegression()
      lr.fit(X, y)
      return RegressionResult(
          predictions=lr.predict(X),
          rsq=lr.score(X,y),
          slope=lr.coef_[1],
          )
    #
    result = self.makeCoStatistic(size)
    title = "Set size: %d, %s" % (size, title)
    col_dict = {
        cn.RATE: cn.RNG, 
        cn.YIELD: cn.RNG,
        cn.MIN: cn.MIN, 
        cn.MAX: cn.MAX, 
        }
    rsqs = []
    slopes = []
    for key, col in col_dict.items():
      if key == cn.MAX:
        new_key = cn.MIN
      else:
        new_key = key
      df = result.get(new_key)
      df = df.sort_values(col)
      plt.scatter(df[col], df[cn.COUNT_MUTATIONS])
      regression_result = regress(df[col], df[cn.COUNT_MUTATIONS])
      rsqs.append(regression_result.rsq)
      slopes.append(regression_result.slope)
      plt.plot(df[col], regression_result.predictions, c='r')
      plt.xlabel(key)
      plt.ylabel("Count")
      new_title = "%s, RSQ=%1.3f, slope=%1.3f" % (
          title, regression_result.rsq, regression_result.slope)
      plt.title(new_title)
      if self._is_plot:
        plt.show()
    df = pd.DataFrame({
        cn.RSQ: rsqs,
        cn.SLOPE: slopes,
        cn.XAXIS: col_dict.keys(),
        })
    df[cn.VALUE] = size
    return df

  @classmethod
  def makeSlopeDF(cls, lines=None, mutation_columns=None, set_sizes=None):
    """
    Constructs dataframes with slopes for lines.
    :return pd.DataFrame:
       cn.LINE, cn.RSQ, cn.SLOPE, cn.VALUE (set size),
       cn.MUTATION_COLUMN,
       cn.XAXIS (cn.RATE, cn.YIELD, cn.MIN, cn.MAX)
    """
    if lines is None:
      lines = [cn.LINE_HA2, cn.LINE_HR2, cn.LINE_UE3]
    if mutation_columns is None:
      mutation_columns = cn.MUTATION_COLUMNS
    if set_sizes is None:
      set_sizes = range(3, 9)
    specification = {
        cn.LINE: lines,
        cn.MUTATION_COLUMN: mutation_columns,
        cn.VALUE: set_sizes,
        }
    dfs = []
    for context in nextStudyContext(specification):
      constraints = [lambda r: r[cn.LINE] == context.line]
      cooccur = cls(context.mutation_column,
          constraints=constraints, is_plot=False)
      df = cooccur.plot(context.value)
      df[cn.LINE] = context.line
      df[cn.MUTATION_COLUMN] = context.mutation_column
      dfs.append(df)
    return pd.concat(dfs)

  @classmethod
  def makeLineCoStatistic(cls, study_context, rc_vector=None):
    """
    Makes statistics for the line and the range constraint vector
    for the range of sizes of isolates present.
    :param StudyContext study_context: specifies
        line, mutation_column
    :param RangeConstraintVector rc_vector:
    :return CoStatistic:
    """
    if study_context.line == cn.LINE_ALL:
        constraints = None
    else:
      constraints = [lambda r: r[cn.LINE] == study_context.line]
    provider = ModelDataDualProvider(study_context.mutation_column,
        constraints=constraints, rc_vector=rc_vector)
    provider.do()
    m_c = cls(mutation_column=study_context.mutation_column,
        provider=provider)
    #
    max_size = len(provider.df_X)  # Number of isolates present
    result = CoStatistic()
    for size in range(2, max_size+1):
      result.concat(m_c.makeCoStatistic(size))
    return result
