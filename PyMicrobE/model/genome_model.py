"""
Model of genome characteristics to predict rate, yield
via regression or classification.
- Acquires data to constructed models.
- Filters observations and does limited transformations
  of predictor variables.
"""


import __init__
import util_data as ud
import constants as cn
from model_data_provider import ModelDataProvider
import util_data as ud
from group_splitter import GroupSplitter
from cv_regression import CVLinearRegression,  \
    CVLassoRegression,  \
    CVForwardRegression, CVBinaryTreeRegression
from cv_classification import CVBinaryTreeClassification,  \
    CVClassification
from isolate_regression import IsolateRegression
from predictor_transformer import PredictorTransformer
from mutation_context import MutationContext
import util

import inspect
import numpy as np
import pandas as pd
from sklearn import linear_model
from sklearn.model_selection import cross_val_predict, cross_val_score
from sklearn.metrics import r2_score

NUM_FOLDS = 31  # Size of the cross validation leave-outs
MAX_STD = 3.0  # Maximum residual STD for selecting outlier cultures
MIN_DEG_FREE = 10  # Minimum degrees of freedom after parameter
                   # estimation and cross validation
NUM_FOLDS_DEFAULT = -1
# Percentile that defines a 0 class; 100-percentile defines a 1 class
PERCENTILE_THRESHOLD = 50.0
MAX_PREDICTOR_COLUMNS = 100  # Maximum number of predictor columns
#


class GenomeModel(object):
  """
  Instances do a single regression for a dependent variable,
  mutation_column (e.g., KEY_MUTATION, POSITION), data constraints.
  Identifies correlated mutations and constructs mutation groups.
  Key methods:
    fit - estimate the parameters of the regression
    score - score the model quality
    predict - provide predicted values
  Note: Must do fit before score and predict
  Currently filters mutations so only select those in multiple lines.
  """

  def __init__(self, depvar, mutation_column, 
    constraints=None,
    transform_type=cn.TRANSFORM_NONE,
    num_folds=NUM_FOLDS_DEFAULT,
    percentile_threshold=PERCENTILE_THRESHOLD,
    is_standardize=True,
    cv_model_cls=CVBinaryTreeClassification,
    **cv_model_flags):
    """
    :param str depvar: dependent variable
    :param str mutation_column:  Name of the column with mutations
    :param list-of-booleanFunction constraints: constraints are predicates
    :param str transform_type: how the predictor variables are transformed
    :param int num_folds: Number of folds in cross validation. If num_folds=-1,
                          then the number of folds is the number of samples
    :param float percentile_threshold: lower threshold values are 0; upper are 1
    :param bool is_standardize: standardize the dependent variable
    :param type(CVModel) cv_model_cls: Model used to analyze data
    :param dict cv_model_flags: Flags passed to the model
    """
    cls = self.__class__
    self._depvar = depvar
    self._mutation_column = mutation_column
    self._constraints = util.setNoneList(constraints)
    self._num_folds = num_folds
    self._cv_model_cls = cv_model_cls
    self._cv_model_flags = cv_model_flags
    self._percentile_threshold = percentile_threshold
    # Get the data
    self.provider = ModelDataProvider(
        MutationContext(self._depvar, self._mutation_column),
        constraints=constraints, is_standardize=is_standardize)
    self.provider.do(transform_type=transform_type)
    self.df_y = self.provider.df_y
    self.df_X = self.provider.df_X
    self.col_y = self.df_y.columns[0]  # Data column
    #
    self.dfs_fit = None  # Dataframe from fit
    self.scores = None  # Vectors of scores
    self.df_predict = None  # Dataframe for prediction

  @classmethod
  def doByLine(cls, lines, depvar, mutation_column, **kwargs):
    """
    Does an analysis by line
    :param list-str lines:
    :return dict: key: line; value: cvRegression
    """
    result = {}
    for line in lines:
      constraint = lambda r: line == r[cn.LINE]
      model = GenomeModel(depvar, mutation_column, 
          constraints=[constraint],
          transform_type=cn.TRANSFORM_LOW_FREQUENCY_ISOLATES,
          **kwargs)
      result[line] = model.fit()
    return result

  def fit(self):
    """
    Constructs the groups used to make cross validation folds and
    runs the cross validation. The resulting fit estimates parameters,
    provides predictions, and scores the regression quality.
    :return CVModel:
       df_parameter, df_predict, df.score
    """
    def cleanDF(df):
      del df[cn.KEY_ISOLATE_DVH]
      del df[cn.KEY_ISOLATE_MMP]
    # Transform data for classification if needed
    if CVClassification in inspect.getmro(self._cv_model_cls):
      ud.makeClassificationData(self.df_X, self.df_y, 
          self.col_y, self._percentile_threshold)
    # Construct the groups
    df_X = self.df_X.reset_index()
    df_y = self.df_y.reset_index()
    df_group = df_X[[cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP]].copy()
    df_group = util.cleanDF(df_group, is_reset_index=True)
    cleanDF(df_X)
    cleanDF(df_y)
    # Check for default action of leave out 1
    if self._num_folds == NUM_FOLDS_DEFAULT:
      num_folds = len(df_y)
    else:
      num_folds = self._num_folds
    # Run the model 
    g_splitter = GroupSplitter(df_X, df_y, df_group, num_folds=num_folds)
    cv_model = self._cv_model_cls(g_splitter, **self._cv_model_flags)
    cv_model.fit()
    return cv_model
