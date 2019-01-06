"""Provides data for models of prediction phenotype from genotype."""


import __init__
import util_data as ud
import constants as cn
from isolate import Isolate
from predictor_transformer import PredictorTransformer
from study_context import StudyContext
import util

import numpy as np
import pandas as pd

PERCENTILE_THRESHOLD = 50.0
MAX_PREDICTOR_COLUMNS = 100  # Maximum number of predictor columns
IS_DEBUG = True


################## CLASSES ###########################
class PredictorProvider(object):
  pass


class ModelDataProvider(PredictorProvider):
  """
  Provides predictor and dependent variables based on constraints
  and specifications of data transformations.
  Data are filtered for outlier cultures. Data are aggregated
  by isolate pair.
  The results are the instance variables df_X, df_y
  """
  #cn.KEY_CULTURE, cn.LINE, cn.GENE_ID, cn.GGENE_ID,
  #cn.POSITION, cn.KEY_MUTATION, cn.RATE, cn.YIELD,
  #cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP
  df_data = None  # Base data

  def __init__(self, context, 
      is_standardize=True, constraints=None, rc_vector=None,
      is_standardize_by_line=False, **kwargs):
    """
    :param MutationContext context: instance variables are:
        depvar, mutation_column
    :param bool is_standardize: standardize the data
    :param list-of-booleanFunction constraints: constraints are predicates on cls.df_data
    :param RangeConstraintVector rc_vector: Range constraint on
        dependent variables
    :param bool is_standardize: data are standardized
    :param bool is_standardize_by_line: data are standardized by line
    """
    cls = self.__class__
    self.context = context
    self._constraints = util.setNoneList(constraints)
    self._rc_vector = rc_vector
    self._is_standardize = is_standardize
    self._is_standardize_by_line = is_standardize_by_line
    #
    if cls.df_data is None:
      cls.df_data = ud.makeIsolateData()
    #
    self.df_X = None
    self.df_y = None
    self.df_y_std = None

  @staticmethod
  def getIsolatesFromIndices(indices):
    """
    Extracts the isolates from the indices of a df_X.
    :param pandas.index indices: 
        cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP
    :return dict: keyed by cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP
      values correspond to rows element in the index
    """
    keys = [n for n in indices.names]
    result = {}
    for idx, key in enumerate(keys):
      result[key] = [v[idx] for v in indices.values]
    return result

  @classmethod
  def getLinesForRows(cls, df):
    """
    Gets the lines for rows in the data.
    :return list-str:
    """
    isolate_dict = cls.getIsolatesFromIndices(df.index)
    return [Isolate.create(v).line for v in
        isolate_dict[cn.KEY_ISOLATE_DVH]]

  def _makeXyDF(self):
    """
    Updates state: self.df_X, self.df_y
    Columns of self.df_X should be strings
    """
    cls = self.__class__
    df = util.selectRows(cls.df_data, self._constraints)
    if len(df) == 0:
      raise ValueError("No data returned with constraints.")
    #
    self.df_y = pd.DataFrame(df[[cn.KEY_CULTURE, self.context.depvar]])
    self.df_y = self.df_y.drop_duplicates()
    self.df_y = self.df_y.set_index(cn.KEY_CULTURE)
    #
    df[cn.COUNT] = 1
    df = util.cleanDF(df, is_reset_index=True)
    self.df_X = util.makeMatrix(df, row_name=cn.KEY_CULTURE,
        column_name=self.context.mutation_column)
    for col in self.df_X.columns:
      self.df_X.rename(columns={col: str(col)}, inplace=True)

  def _aggregateByIsolatePair(self):
    """
    Creates X, y aggregated by and indexed by isolate pair.
    :param pd.DataFrame df: indexed by culture
    :return pd.DataFrame: indexed by cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP
    """
    cls = self.__class__
    # Map from culture to isolate pair
    df_merge = cls.df_data[
        [cn.KEY_CULTURE, cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP]]
    df_merge = df_merge.copy()
    df_merge = util.cleanDF(df_merge, is_reset_index=True)
    #
    def aggregate(df, aggregate_type=cn.AVG):
      df_result = df.copy()
      df_result = df_result.reset_index()
      df_result = df_result.merge(df_merge, on=cn.KEY_CULTURE, how='inner')
      del df_result[cn.KEY_CULTURE]
      df_grp  = df_result.groupby([cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP])
      if aggregate_type == cn.AVG:
        df_result = df_grp.mean()
      elif aggregate_type == cn.STD:
        df_count = df_grp.count()
        df_count = df_count.applymap(lambda v: np.sqrt(v))
        df_result = df_grp.std()
        df_result = df_result.divide(df_count, fill_value=0)
      else:
        raise ValueError("Invalid aggreate_type %s" % aggregate_type)
      df_result = df_result.copy()
      return df_result
    #
    self.df_X = aggregate(self.df_X)
    df_y = self.df_y.copy()
    self.df_y = aggregate(df_y)
    self.df_y_std = aggregate(df_y, aggregate_type=cn.STD)

  def do(self, transform_type=cn.TRANSFORM_NONE, 
      **kwargs):
    """
    Creates the predictor and dependent variable for the regression.
    The columns are:
      X - mutations, line 
      y - depvar
    Both DFs are indexed by cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP.
    :param str transform_type: how the predictor variables are transformed
    :raises ValueError: if empty dataframe with constraints
    Results are in self.df_X, self.df_y
    """
    self._makeXyDF()
    #
    transformer = PredictorTransformer(self.df_X, 
        self.context.mutation_column)
    if transform_type == cn.TRANSFORM_DEFAULT:
      transformer.transformDefault(**kwargs)
    elif transform_type == cn.TRANSFORM_LOW_FREQUENCY_ISOLATES:
      transformer.transformLowFrequencyColumns(MAX_PREDICTOR_COLUMNS)
    elif transform_type == cn.TRANSFORM_ONLY_LOW_FREQUENCY_ISOLATES:
      transformer.transformOnlyLowFrequencyColumns(
          MAX_PREDICTOR_COLUMNS)
    else:
      pass  # no transformation of the predictors
    self.df_X = transformer.df_X
    # Aggregate data by isolate pairs
    self._aggregateByIsolatePair()
    # Standardize the values
    if self._is_standardize_by_line:
      self.df_y = self._standardizeByLine()
    elif self._is_standardize:
      util.standardize(self.df_y)
    # Prune based on aggregated values of dependent variable
    if self._rc_vector is not None:
      indices = self._rc_vector.findSatisfiedRows(self.df_y)
      self.df_X = self.df_X.loc[indices]
      self.df_y = self.df_y.loc[indices]
      self.df_y_std = self.df_y_std.loc[indices]

  def getMutations(self):
    return [str(c) for c in self.df_X.columns]

  def _standardizeByLine(self):
    """
    Standardizes values separately for each line
    :return list-str: lines
    """
    lines = [Isolate.create(v[0]).line 
        for v in self.df_X.index.tolist()]
    df_old = self.df_y.copy()
    df_old[cn.LINE] = lines
    dfs = []
    for line in set(lines):
      df = df_old[df_old[cn.LINE] == line].copy()
      del df[cn.LINE]
      util.standardize(df)
      dfs.append(df)
    return pd.concat(dfs, sort=True)
      

class ModelDataDualProvider(PredictorProvider):
  """Data Provider with both rate and yield."""

  def __init__(self, mutation_column, **kwargs):
    """
    :param MutationContext context:
    :param list-of-booleanFunction constraints: constraints are predicates on cls.df_data
    """
    self.mutation_column = mutation_column
    self._kwargs = kwargs
    self.df_X = None
    self.df_ys = {cn.RATE: None, cn.YIELD: None}
    self.df_y_stds = {cn.RATE: None, cn.YIELD: None}

  def do(self, **kwargs):
    """
    Constructs the providers. Change the name of the df_y column to cn.VALUE.
    """
    df_Xs = {}
    for depvar in cn.DEPVARS:
      provider = ModelDataProvider(
          StudyContext(depvar=depvar, 
          mutation_column=self.mutation_column),
          **self._kwargs)
      provider.do(**kwargs)
      col_y = provider.df_y.columns[0]
      provider.df_y.rename(columns={col_y: cn.VALUE}, inplace=True)
      df_Xs[depvar] = provider.df_X
      self.df_ys[depvar] = provider.df_y
      #
      provider.df_y_std.rename(
          columns={col_y: cn.VALUE}, inplace=True)
      self.df_y_stds[depvar] = provider.df_y_std
    # Reconsile constraints on cn.RATE, cn.YIELD
    index = df_Xs[cn.RATE].index.intersection(df_Xs[cn.YIELD].index)
    self.df_X = df_Xs[cn.RATE].loc[index]
    for key in cn.DEPVARS:
      self.df_ys[key] = self.df_ys[key].loc[index]
      self.df_y_stds[key] = self.df_y_stds[key].loc[index]


################ FUNCTIONS ###################
def makeTransformedData(data_cls=ModelDataDualProvider,
    max_columns=100, **kwargs):
  """
  Makes a data provider with standard transformations.
  :param type data_cls: provider class to create
  :param dict kwargs: arguments passed to create data provider
  :return ModelDataDualProvider:
  """
  def getArgsValueAndDeleteKey(key):
    if key in kwargs:
      value = kwargs[key]
      del kwargs[key]
      return value
    else:
      return None
  #
  mutation_column=getArgsValueAndDeleteKey("mutation_column")
  if mutation_column is None:
    raise ValueError("Must specify mutation_column.")
  context = StudyContext(
      mutation_column=mutation_column,
      depvar=getArgsValueAndDeleteKey("depvar"),
      )
  #
  if data_cls == ModelDataProvider:
    provider = data_cls(context, **kwargs)
  else:
    provider = data_cls(mutation_column, **kwargs)
  provider.do()
  # Filter the mutation data
  transformer = PredictorTransformer(provider.df_X,
      mutation_column)
  transformer.makeMutationGroups()
  transformer.filterLowFrequencyColumns(max_columns)
  provider.df_X = transformer.df_X
  return provider
