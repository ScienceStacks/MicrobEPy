"""
Performs cross validation for an IsolateModel

A model is an instance of IsolateModel.
CrossValidation has the following responsibilities:
  1. Construct estimates and their standard devitations
     for all instances of the IsolateModel class.
  2. Construct predictions and their residuals for all
     model instances of the IsolateModel class.

Cross validation is done in two ways. In the first case,
we leave out a culture that is being predicted.
The second case is to leave out the isolate pair being predicted.
"""

import __init__
import constants as cn
from isolate_model import IsolateModel

import pandas as pd
import numpy as np

KEY_NAMES = [cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP]
MAX_REPLICATION = 4


##########################################
# CLASSES
##########################################
class CVIsolateModel(object):
  """
  Does cross validation for a IsolateModel.
    1. Selects data for IsolateModel
    2. Aggregates results of estimates
    3. Computes standard deviations of estimates
  """

  def __init__(self, model_cls, leave_out=cn.KEY_CULTURE,
      isPermittedRow=lambda r: True, **kwargs):
    """
    :param IsolateModel model_class: A class inheriting from IsolateModel
    :param str leave_out: Indicates what is being left out in the cross validation.
                          cn.KEY_CULTURE (culture) or cn.KEY_ISOLATE (isolate)
    :param function isPermittedRow: row is used if True
        boolean valued function of a dataframe row
    :param dict kwargs: keyword arguments used when instantiating IsolateModel
    """
    if not issubclass(model_cls, IsolateModel):
      raise ValueError("model_cls must inherit from IsolateModel!")
    self._model_cls = model_cls
    self._kwargs = kwargs
    self._leave_out = leave_out
    self._key_values = model_cls.getKeyValues(isPermittedRow, **kwargs)
    self._isPermittedRow = isPermittedRow

  def estimate(self):
    """
    Estimates the parameters of the model and evaluates the model
    Does cross validation for the regression model.
    :return dict:
      cn.AVG: DataFrame of average values
        cn.KEY_ISOLATE_DVH
        cn.KEY_ISOLATE_MMP
        cn.ESTIMATE
        columns: model parameters
      cn.STD: DataFrame of standard deviations
        cn.KEY_ISOLATE_DVH
        cn.KEY_ISOLATE_MMP
        cn.ESTIMATE
        columns: model parameters
      cn.RESIDUAL: DataFrame related to estimates
        cn.KEY_CULTURE
        cn.OBSERVED
        cn.ESTIMATE
        cn.RESIDUAL
        cn.RESIDUALSTD - residuals in units of standard deviation
      cn.RSQ: float - R2 for the model
    """
    #
    def updateEstimatePredict(key_value, dfs_estimate, dfs_predict, idx=None):
      """
      Computes and augments the estimates and predictions.
      Handles case where an index is left out and when it is not left out.
      :param tuple-of-str key_value: DVH, MMP
      :param list-of-dataframe dfs_estimate:
      :param list-of-dataframe dfs_predict:
      :param int idx: integer index of replication of culture
      """
      model = self._model_cls(key_value, 
          isPermittedRow=self._isPermittedRow, 
          leave_out=self._leave_out, **self._kwargs)
      if self._leave_out == cn.KEY_CULTURE:
        isIndex_estimate = lambda x: True
        isIndex_predict = lambda x: True
      else:
        max_idx = model.getDataSize() - 1
        if (idx > max_idx) or (max_idx <= 0):
          return
        isIndex_estimate = lambda x: x % max_idx != idx
        isIndex_predict = lambda x: x % max_idx == idx
      df = model.estimate(isIndex=isIndex_estimate)
      if len(df) > 0:
        dfs_estimate.append(df)
        df_predict = model.predict(isIndex=isIndex_predict)
        if df_predict is not None:  # None if no estimate can be produced
          df_predict[cn.RESIDUAL] = df_predict[cn.OBSERVED] - df_predict[cn.ESTIMATE]
          dfs_predict.append(df_predict)
    #
    result = {}
    dfs_estimate = []
    dfs_predict = []
    # Construct the models and accumulate results
    for key_value in self._key_values:
      if self._leave_out == cn.KEY_CULTURE:
        model = updateEstimatePredict(key_value, dfs_estimate, dfs_predict)
      else:
        for idx in range(MAX_REPLICATION):
          updateEstimatePredict(key_value, dfs_estimate, dfs_predict, idx=idx)
    # Compute statistics
    df_full_estimate = pd.concat(dfs_estimate)
    groupby = df_full_estimate.groupby(KEY_NAMES)
    result[cn.AVG] = groupby.mean().reset_index()
    result[cn.STD] = groupby.std().reset_index()
    df_full_predict = pd.concat(dfs_predict)
    groupby = df_full_predict.groupby(KEY_NAMES).var()
    df = pd.concat(dfs_predict)
    rsq = 1 - np.var(df[cn.RESIDUAL]) / np.var(df[cn.OBSERVED])
    result[cn.RSQ] = min(max(rsq, 0.0), 1.0)
    residual_columns = [cn.OBSERVED, cn.ESTIMATE, cn.RESIDUAL]
    other_columns = set(
        [cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP, cn.KEY_CULTURE]).intersection(df.columns)
    residual_columns.extend(other_columns)
    result[cn.RESIDUAL] = df[residual_columns].copy()
    std = np.std(result[cn.RESIDUAL][cn.RESIDUAL])
    result[cn.RESIDUAL][cn.RESIDUALSTD] = result[cn.RESIDUAL][cn.RESIDUAL] / std
    return result
