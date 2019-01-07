""" Cross validated classification models. """


import microbepy_init
import constants as cn
import binary_tree_classification as btc
from cv_model import CVModel
import util

import copy
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score
from sklearn import linear_model


class CVClassification(CVModel):
  """ Codes common to classification. """
  
  def __init__(self, model, splitter, **model_args):
    """
    :param dict model_args:
                                    a class value.
    """
    super(CVClassification, self).__init__(model, splitter)

  def calcPredictScore(self):
    """
    Computes accuracy from the predict dataframe.
    """
    #
    counts = [1 if r == 0 else 0 
             for r in self.df_predict[cn.RESIDUAL]]
    return sum(counts)/(1.0*len(counts))


################################################################
#              Classification Classes                              #
################################################################

class CVBinaryTreeClassification(CVClassification):
  """
  Methods for ordinary least squares regression
  """

  def __init__(self, splitter, **model_args):
    """
    :param dict model_args: arguments passed in model construction
    """
    model = btc.BinaryTreeClassification(**model_args)
    super(self.__class__, self).__init__(model, splitter)

  def fitDF(self, df_X, df_y):
    """
    Fits single model.
    :param pd.DataFrame df_X:
    :param pd.DataFrame df_y:
    :return pd.DataFrame, model: dataframe has parameter values
                                 for this fit. A 1 means
                                 that a parameter was present.
    """
    fitted_model = self._model.fit(df_X, df_y)
    if fitted_model is None:
      import pdb; pdb.set_trace()
    df_parameter = fitted_model.getParameterDF()
    missing_columns = set(df_X.columns).difference(
        df_parameter.columns)
    for col in missing_columns:
      df_parameter[col] = 0
    df_parameter[cn.ACC] = fitted_model.score(df_X, df_y)
    return df_parameter, fitted_model

  def predictDF(self, df_X, df_y):
    """
    DataFrame of predictions across cross validation instances.
    :param pd.DataFrame df_X:
    :param pd.DataFrame df_y:
    :return pd.DataFrame: cn.ESTIMATE, cn.OBSERVED, cn.RESIDUAL
    """
    cls = self.__class__
    super(cls, self)._predictDF(df_X, df_y)
    #
    return self._fitted_model.predictDF(df_X, df_y)
