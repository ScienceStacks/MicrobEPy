"""Abstract class for cross validation for models, either regression or classification."""

from microbepy.common import constants as cn
from microbepy.common import util

import copy
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score


class CVModel(object):
  """
  Abstract class for cross validation.
  Results are:
    self.df_parameter - cn.AVG, cn.STD, cn.COUNT, indexed by parameter name
    self.df_predict - cn.ESTIMATE, cn.OBSERVED, cn.RESIDUAL
    self.score - value in [0, 1]; 0 is low; 1 is high.
  """

  def __init__(self, model, splitter):
    """
    :param sklearn.linear_model model:
    :param Iterator splitter: 
    """
    self._model = model  #  Instantiated but unfitted model
    self._splitter = splitter
    # Outputs
    # State produced by fit method
    self.score = None  # CV score
    # Parameter estimates. Indexed by parameter name,
    # cn.AVG, cn.STD
    self.df_parameter = None  # Parameter dataframe
    # Indexed as predictor variables
    # cn.ESTIMATE, cn.OBSERVED, cn.RESIDUAL
    self.df_predict = None  # Parameter dataframe
    # Results of individual models
    self.df_parameter_accum = None
    # Models constructed
    self._fitted_model = None  # Model after a parameter fitting
    self.fitted_models = []  # All models used

  ################################################################
  #              Methods Overriden by Subclasses                 #
  ################################################################
  def fitDF(self, df_X, df_y):
    """
    Fit a model.
    :param pd.DataFrame df_X:
    :param pd.DataFrame df_y:
    :return pd.DataFrame, model: parameter estimates, fitted model
    """
    raise RuntimeError("Must override.")

  def predictDF(self, df_X, df_y):
    """
    Performs prediction for an existing fit.
    :param pd.DataFrame df_X:
    :param pd.DataFrame df_y:
    :return pd.DataFrame: parameter estimates
    This is a default method that may be overridden.
    """
    cls = self.__class__
    self._predictDF(df_X, df_y)
    try:
      y_values = [v[0] for v in self._fitted_model.predict(df_X)]
    except:
      import pdb; pdb.set_trace()
    df = cls._makePredictDF(df_y, y_values)
    return df

  ################################################################
  #              Exposed Methods                                 #
  ################################################################
  def fit(self, num_residual_leaveouts=0):
    """
    Fits the model. Produces the score, parameters, and predictions.
    State affected: sef.score, self.df_parameter, self.df_predict
    :param int num_residual_leaveouts: Number of extreme residuals to leaveout
    """
    self.df_parameter_accum = pd.DataFrame()  # Accumulate parameter estimates
    self.df_predict = pd.DataFrame()  # Predictions
    num_folds = 0
    for dfs in self._splitter:
      num_folds += 1  # Counts the number of folds used in the calculation
      df_parameter, self._fitted_model = self.fitDF(
          dfs[cn.TRAIN_X], dfs[cn.TRAIN_Y])
      self.fitted_models.append(self._fitted_model)
      columns = df_parameter.columns.tolist()
      # Only predict if there's a non-null result
      if len(columns) > 1:  
        for col in [cn.RSQ, cn.ACC]:
          if col in columns:
            columns.remove(col)
        self.df_parameter_accum = util.appendWithColumnUnion(
            self.df_parameter_accum, df_parameter)
        #
        df_predict = self.predictDF(dfs[cn.TEST_X], dfs[cn.TEST_Y])
        self.df_predict = self.df_predict.append(df_predict)
    df_mean = self.df_parameter_accum.mean()
    df_std = self.df_parameter_accum.std()
    self.df_parameter = pd.DataFrame({
        cn.AVG: df_mean,
        cn.STD: df_std,
        cn.COUNT: np.repeat(num_folds, len(df_mean)),
        }) 
    # Adjust for residual leaveouts
    if num_residual_leaveouts > 0:
      df = self.df_predict.copy()
      df[cn.RESIDUAL] = df[cn.RESIDUAL].apply(lambda v: np.abs(v))
      df.sort_values(cn.RESIDUAL, ascending=False)
      cultures = df.index[0:num_residual_leaveouts]
      rows = [r for c,r in self.df_predict.iterrows() 
              if not c in cultures]
      self.df_predict = pd.DataFrame(rows)
    # Score the cross validation
    if len(self.df_predict) == 0:
      self.score = 0.0
    else:
      self.score = self.calcPredictScore()

  ################################################################
  #              Methods used by subclasses                      #
  ################################################################
  @staticmethod
  def _getParameterDF(df_X, values):
    """
    Constructs a DataFrame with values corresponding to columns.
    :param pd.DataFrame df:X:
    :param list-float values:
    :return pd.DataFrame: same columns as self.df_X
    """
    rows = {}
    for n,col in enumerate(df_X.columns):
      rows[col] = [values[n]]
    df = pd.DataFrame(rows)
    return df

  def _predictDF(self, df_X, df_y):
    """
    Performs prediction for an existing fit.
    :param pd.DataFrame df_X:
    :param pd.DataFrame df_y:
    :return pd.DataFrame: parameter estimates
    """
    if self._fitted_model is None:
      raise ValueError("Must do fit before prediction!")

  @staticmethod
  def _makePredictDF(df_y, y_values):
    """
    Constructs prediction dataframe.
    :param pd.DataFrame df_y:
    :param list-values y_values:
    :return pd.DataFrame: cn.ESTIMATE, cn.OBSERVED, cn.RESIDUAL
    """
    column = df_y.columns[0]
    df = pd.DataFrame({cn.OBSERVED: df_y[column]})
    df[cn.ESTIMATE] = y_values
    df[cn.RESIDUAL] = df[cn.OBSERVED] - df[cn.ESTIMATE]
    return df

  def findSignificantParameters(self, 
      predicate=lambda r: not np.isclose(abs(r[cn.AVG]), 0)):
    """
    Finds the parameters of the model that are significant
    according to a predicate.
    :param BooleanFunction predicate: function of a df_parameter row
    :return list-str: parameter names
    """
    parameters = [i for i, r in self.df_parameter.iterrows()
                  if predicate(r)]
    for name in [cn.RSQ, cn.ACC]:
      if name in parameters:
        parameters.remove(name)
    return parameters
