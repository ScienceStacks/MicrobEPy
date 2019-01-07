""" Cross validated regression models. """


import microbepy_init
import constants as cn
import binary_tree_regression as btr
from cv_model import CVModel
import util

import copy
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score
from sklearn import linear_model


class CVRegression(CVModel):
  """ Codes common to regressions. """

  def calcScore(self, df_X, df_y):
    """
    Constructs the score for an existing fit.
    :param pd.DataFrame df_X:
    :param pd.DataFrame df_y:
    :return float:
    """
    rsq = self._model.score(df_X, df_y)
    return rsq

  def calcPredictScore(self):
    """
    Computes RSQ from the predict dataframe.
    """
    return r2_score(self.df_predict[cn.OBSERVED].tolist(),
          self.df_predict[cn.ESTIMATE].tolist())


################################################################
#              Regression Classes                              #
################################################################

class CVLinearRegression(CVRegression):
  """
  Methods for ordinary least squares regression
  """

  def __init__(self, splitter, **model_args):
    """
    :param dict model_args: arguments passed in model construction
    """
    model = linear_model.LinearRegression(**model_args)
    super(CVRegression, self).__init__(model, splitter)

  def fitDF(self, df_X, df_y):
    """
    Fit a Linear regression model.
    :param pd.DataFrame df_X:
    :param pd.DataFrame df_y:
    :return pd.DataFrame, model:
    """
    cls = self.__class__
    if not util.isColinear(df_X):
      fitted_model = self._model.fit(df_X, df_y)
      values = fitted_model.coef_[0]
      df = cls._getParameterDF(df_X, values)
    else:
      fitted_model = None
      print ("Colinear fold")
      df = pd.DataFrame()
    df[cn.RSQ] = self.calcScore(df_X, df_y)
    return df, fitted_model


class CVLassoRegression(CVRegression):
  """
  Methods for Lasso regression
  """

  def __init__(self, splitter, **model_args):
    """
    :param dict model_args: arguments passed in model construction
    """
    model = linear_model.Lasso(**model_args)
    super(CVRegression, self).__init__(model, splitter)

  def fitDF(self, df_X, df_y):
    """
    Fit a Linear regression model.
    :param pd.DataFrame df_X:
    :param pd.DataFrame df_y:
    :return pd.DataFrame, model:
    """
    cls = self.__class__
    if not util.isColinear(df_X):
      fitted_model = self._model.fit(df_X, df_y)
      df = cls._getParameterDF(df_X, self._model.coef_.tolist())
    else:
      self._fitted_model = None
      print ("Colinear fold")
      df = pd.DataFrame()
    df[cn.RSQ] = self.calcScore(df_X, df_y)
    return df, fitted_model

  def predictDF(self, df_X, df_y):
    """
    Performs prediction for an existing fit.
    :param pd.DataFrame df_X:
    :param pd.DataFrame df_y:
    :return pd.DataFrame: parameter estimates
    """
    cls = self.__class__
    super(cls, self)._predictDF(df_X, df_y)
    #
    y_values = self._fitted_model.predict(df_X).tolist()
    df = cls._makePredictDF(df_y, y_values)
    return df


class CVForwardRegression(CVRegression):
  """
  Implements forward stagewise regression.
  Use a LinearRegression model in the constructor.
  """
  
  def __init__(self, splitter, **model_args):
    model = linear_model.LinearRegression(**model_args)
    super(CVRegression, self).__init__(model, splitter)
    self.selected_columns = []  # Columns selected for the model

  def fitDF(self, df_X, df_y):
    """
    Fit a Linear regression model.
    :param pd.DataFrame df_X:
    :param pd.DataFrame df_y:
    :return pd.DataFrame: parameter estimates
    """
    cls = self.__class__
    MIN_CHANGE = 0.1  # Minimal increase in R2
    MIN_DEG_FREE = 5
    score = 0.0
    remaining_columns = df_X.columns.tolist()
    depvar = df_y.columns[0]
    done = False
    df_X_cur = pd.DataFrame()
    df_res = df_y.copy()
    fitted_model = self._model
    while not done and (len(remaining_columns) > 0):
      if len(df_y) - len(df_X_cur.columns) < MIN_DEG_FREE:
        break
      temp_model = copy.deepcopy(fitted_model)
      df = df_X[remaining_columns]
      df_corr = util.correlateWithPredictors(df, df_res)
      column = df_corr.index[0]
      df_X_cur[column] = df_X[column]
      temp_model = temp_model.fit(df_X_cur, df_y)
      new_score = temp_model.score(df_X_cur, df_y)
      if (len(df_X_cur.columns) > 1)  \
           and util.isColinear(df_X_cur):
        del df_X_cur[column]
        remaining_columns.remove(column)
        continue
      if new_score < score + MIN_CHANGE:
        del df_X_cur[column]
        break
      fitted_model = copy.deepcopy(temp_model)
      score = new_score
      remaining_columns.remove(column)
      predicteds = fitted_model.predict(df_X_cur)
      df_res[depvar] = df_y[depvar] - [v[0] for v in predicteds]
    if len(df_X_cur.columns) > 0:
      values = fitted_model.coef_[0]
      df_parameter = cls._getParameterDF(df_X_cur, values)
      self.selected_columns = df_X_cur.columns.tolist()
    else:
      df_parameter = pd.DataFrame()
      self.selected_columns = []
    df_parameter[cn.RSQ] = score
    return df_parameter, fitted_model

  def predictDF(self, df_X, df_y):
    """
    Performs prediction for an existing fit.
    :param pd.DataFrame df_X:
    :param pd.DataFrame df_y:
    :return pd.DataFrame: parameter estimates
    """
    cls = self.__class__
    super(cls, self)._predictDF(df_X, df_y)
    #
    df_XX = df_X[self.selected_columns]
    try:
      y_values = [v[0] for v in self._fitted_model.predict(df_XX)]
    except:
      import pdb; pdb.set_trace()
    df = cls._makePredictDF(df_y, y_values)
    return df


class CVBinaryTreeRegression(CVRegression):

  def __init__(self, splitter, **model_args):
    model = btr.BinaryTreeRegression(**model_args)
    super(self.__class__, self).__init__(model, splitter)

  def fitDF(self, df_X, df_y):
    """
    :param pd.DataFrame df_X:
    :param pd.DataFrame df_y:
    :return pd.DataFrame, model: parameter estimates, cn.RSQ
      Parameter for each column. Value is 1 if present
      in the Tree; 0 otherwise.
    """
    fitted_model = copy.deepcopy(self._model)
    fitted_model.fit(df_X, df_y)
    nodes = fitted_model.findNodes()
    split_columns = [n.split_column for n in nodes 
                    if n.split_column != cn.LEAF_NODE]
    rows = {}
    for col in df_X.columns:
      if col in split_columns:
        rows[col] = [1]
      else:
        rows[col] = [0]
    df = pd.DataFrame(rows)
    df[cn.RSQ] = fitted_model.score(df_X, df_y)
    #
    return df, fitted_model

  def predictDF(self, df_X, df_y):
    """
    :param pd.DataFrame df_X:
    :param pd.DataFrame df_y:
    :return pd.DataFrame: cn.ESTIMATE, cn.OBSERVED, cn.RESIDUAL
    """
    cls = self.__class__
    super(cls, self)._predictDF(df_X, df_y)
    #
    return self._fitted_model.predictDF(df_X, df_y)
