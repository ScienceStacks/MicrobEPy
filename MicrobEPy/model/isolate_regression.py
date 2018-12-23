"""
Performs a regression for a pair of isolates, one a DVH and the other
an MMP. The DVH isolate is indexed by i; the MMP isolate is
indexed by j.

The model is that there is an effect due to DVH j, denoted by
\alpha_i, an effect due to MMP j, denoted by \gamma_j, and
an effect due to their interaction, denoted by \gamma_ij.

The model is:
  y_ijk = alpha_i + gamma_j + gamma_ij + epsilon_ijk,
where k indexes the culture.

Thus, in all cases the dependent variables are the
of DVH and MMP isolates.

The different models estimate alpha_i, gamma_j, gamma_ij in
different ways.

Data are identified at two levels:
 - Observeations are identified by cn.KEY_CULTURE
 - Cases are isolate pairs and are identified by
   cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP
"""

import __init__
import constants as cn
import isolate_model as im
from cv_isolate_model import CVIsolateModel
from isolate import Isolate
import util

import collections
from sklearn import linear_model
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

COEF_DVH = 'coef_dvh'
COEF_MMP = 'coef_mmp'


##########################################
# FUNCTIONS
##########################################
def makeAvgAndStdDF(groupby):
  """
  Calculates the mean and its standard deviation
  :param pd.DataFrameGroupBy groupby:
  :return pd.DataFrame: INDEX, cn.AVG, cn.STD
  """
  counts = groupby.count().apply(np.sqrt).values.flatten()
  pairs = zip(counts, groupby.std().values.flatten())
  stds = [v/c for c,v in pairs]
  df = pd.DataFrame({
      cn.KEY_ISOLATE: list(groupby.indices.keys()),
      cn.AVG: groupby.mean().values.flatten(), 
      cn.STD: stds,
      })
  df.set_index(cn.KEY_ISOLATE, inplace=True)
  return df


##########################################
# CLASSES
##########################################
class IsolateRegression(im.IsolateModel):
  """
  Abstract class with code common to isolate regressions.
  """

  def __init__(self, isolate_pair, depvar=cn.RATE, 
      isPermittedRow=im.ISPERMITTEDROW, **kwargs):
    """
    :param str depvar: Dependent variable for regression
                        cn.RATE or cn.YIELD
    :param list-of-str isolate_pair: list of isolate string
    :param function isPermittedRow: boolean function of rows
      of df_coculture
    """
    cls = self.__class__
    self._isolate_mmp = [iso for iso in isolate_pair
                         if cn.SPECIES_MIX_MMP in iso][0]
    self._isolate_dvh = [iso for iso in isolate_pair
                         if cn.SPECIES_MIX_DVH in iso][0]
    self._depvar = depvar
    self._isPermittedRow = isPermittedRow
    self.df_estimate = None
    self.df_predict = None
    cls._getData()

  @classmethod
  def getKeyValues(cls, isPermittedRow=im.ISPERMITTEDROW, **kwargs):
    """
    :return list-of-list-of-str: cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP
    """
    cls._getData()
    df = cls.dfs_coculture[cn.RATE]
    sel = [isPermittedRow(r) for _,r in df.iterrows()]
    df_result = df.loc[sel, [cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP]]
    df_result = df_result.drop_duplicates()
    result = [[r[cn.KEY_ISOLATE_DVH], r[cn.KEY_ISOLATE_MMP]] 
              for _,r in df_result.iterrows()]
    return result

  def _getDepvars(self, isIndex=lambda v: True):
    """
    Obtains values of the dependent variable for the given DVH, MMP
    :param Function isIndex: True of data index is to be included
    :return pd.DataFrame: cn.DEPVAR, cn.KEY_CULTURE
    """
    cls = self.__class__
    constraint = lambda r: (
        (r[cn.KEY_ISOLATE_DVH] == self._isolate_dvh) and
        (r[cn.KEY_ISOLATE_MMP] == self._isolate_mmp) and
        self._isPermittedRow(r)
        )
    df = cls.getDepvarDF(self._depvar, isIndex, constraint=constraint)
    return df

  def getDataSize(self):
    """
    :return int: Number of replications
    """
    cls = self.__class__
    df = cls.dfs_coculture[self._depvar]
    count = sum([1 if
        (r[cn.KEY_ISOLATE_DVH] == self._isolate_dvh) and
        (r[cn.KEY_ISOLATE_MMP] == self._isolate_mmp) 
        else 0 
        for _,r in df.iterrows()
        if self._isPermittedRow(r)])
    return count

  @classmethod
  def getDepvarDF(cls, depvar, isIndex=lambda v: True, 
      constraint=lambda v: True):
    """
    Obtains values of the dependent variable for the constraints
    :param str depvar: dependent variable
    :param Function isIndex: True of data index is to be included
    :param Function constraint: boolean function of df_coculture row
    :return pd.DataFrame: cn.DEPVAR, cn.KEY_CULTURE
    """
    df_coculture = cls.dfs_coculture[depvar]
    sel = [constraint(r) for _,r in df_coculture.iterrows()]
    df = df_coculture.loc[sel, :]
    depvars = df[cn.DEPVAR].tolist()
    cultures = df[cn.KEY_CULTURE].tolist()
    rows = {cn.DEPVAR: [], cn.KEY_CULTURE: []}
    for idx in range(len(depvars)):
      if isIndex(idx):
        rows[cn.DEPVAR].append(depvars[idx])
        rows[cn.KEY_CULTURE].append(cultures[idx])
    return pd.DataFrame(rows)

  def predict(self, isIndex=lambda x: x == 0):
    """
    Constructs a DataFrame for the predictions
    :param Function isIndex:
    :return pd.DataFrame:
      cn.KEY_CULTURE
      cn.KEY_ISOLATE_DVH
      cn.KEY_ISOLATE_MMP
      cn.ESTIMATE - estimated value
      cn.OBSERVED
    """
    if self.df_estimate is None:
      raise ValueError("Must do estimate before doing predict.")
    if len(self.df_estimate.index) == 0:
      return None
    df = self._getDepvars(isIndex=isIndex)
    if len(df) == 0:
        return None
    self.df_predict = pd.DataFrame({
        cn.OBSERVED: df[cn.DEPVAR],
        cn.KEY_CULTURE: df[cn.KEY_CULTURE],
        })
    self.df_predict[cn.ESTIMATE] = self.df_estimate[cn.ESTIMATE].tolist()[0]
    self.df_predict[cn.KEY_ISOLATE_DVH] = self._isolate_dvh
    self.df_predict[cn.KEY_ISOLATE_MMP] = self._isolate_mmp
    return self.df_predict

  @classmethod
  def makeEstimateDFS(cls, cvsize=3, depvar=cn.RATE,
      isPermittedRow=im.ISPERMITTEDROW, **kwargs):
    """
    :param int cvsize: size of the cross validation set
    :param str depvar: dependent variable
    :param function isPermittedRow: True if row is permitted
        allows for filtering
    :return dict: Cross validation dataframes
    """
    cross = CVIsolateModel(cls, depvar=depvar,
        isPermittedRow=isPermittedRow, **kwargs)
    result = cross.estimate()
    return result

  @classmethod
  def getSmallResidualCultures(cls, max_std=3, **kwargs):
    """
    Finds the cultures for those points with small absolute value of residuals.
    :param arguments kwargs: arguments for makeEstimateDFS
    :param int max_std: Maximum absolute value of the residual in units of standard deviation
    :return list-of-str: list of culures
    """
    dfs = NonParametricIsolateRegression.makeEstimateDFS(kwargs)
    df = dfs[cn.RESIDUAL]
    cultures = [r[cn.KEY_CULTURE] for _,r in df.iterrows() 
                if (r[cn.RESIDUALSTD] < max_std) 
                and (r[cn.RESIDUALSTD] > -max_std)]
    return cultures


class ParametricIsolateRegression(IsolateRegression):

  @classmethod
  def getKeyValues(cls, isPermittedRow=im.ISPERMITTEDROW, **kwargs):
    cls._getData()
    parent_cls = util.getParentClass(cls)
    initial_isolates = parent_cls.getKeyValues(
        isPermittedRow=isPermittedRow)
    alpha_isolates = set(cls.dfs_ancestral[cn.RATE].index)
    result =  [ (d, m) for (d, m) in initial_isolates
             if (d in alpha_isolates) and (m in alpha_isolates)]
    if len(result) == 0:
      import pdb; pdb.set_trace()
    return result

  def estimate(self, isIndex=lambda x: True, **kwargs):
    """
    Estimates the mean and standard devition of values of the model constants.
    :param Function isIndex:
    :return pd.DataFrame:
      cn.KEY_ISOLATE_DVH
      cn.DEPVAR
      im.ALPHA_AVG
      im.ALPHA_STD
      im.BETA_AVG
      im.BETA_STD
      im.GAMMA_AVG
      im.GAMMA_STD
      cn.COUNT - number replications of the isolate pair
      cn.ESTIMATE - estimated value for the combination of evolved isolates
    """
    cls = self.__class__
    cls._getData()
    #
    df_ancestral = cls.dfs_ancestral[self._depvar]
    alpha_avg = df_ancestral.loc[self._isolate_dvh, cn.AVG]
    alpha_std = df_ancestral.loc[self._isolate_dvh, cn.STD]
    beta_avg = df_ancestral.loc[self._isolate_mmp, cn.AVG]
    beta_std = df_ancestral.loc[self._isolate_mmp, cn.STD]
    #
    depvars = self._getDepvars(isIndex=isIndex)[cn.DEPVAR]
    gamma_avg = np.mean(depvars) - alpha_avg - beta_avg
    # TODO: Calculate least squares estimator of STD
    gamma_std = np.std(depvars)
    estimate = alpha_avg + beta_avg + gamma_avg
    self.df_estimate = pd.DataFrame({
        cn.KEY_ISOLATE_DVH: [self._isolate_dvh],
        cn.KEY_ISOLATE_MMP: [self._isolate_mmp],
        cn.DEPVAR: [self._depvar],
        im.ALPHA_AVG: [alpha_avg],
        im.ALPHA_STD: [alpha_std],
        im.BETA_AVG: [beta_avg],
        im.BETA_STD: [beta_std],
        im.GAMMA_AVG: [gamma_avg],
        im.GAMMA_STD: [gamma_std],
        cn.COUNT: [self.getDataSize()],
        cn.ESTIMATE: [estimate],
        })
    return self.df_estimate

  @classmethod
  def plotAlphaEstimate(cls, depvar, is_test=False, plot_parms=None):
    """
    Scatter plot with DVH, MMP alphas and sizes of estimates.
    :param str depvar: dependent variable (cn.RATE, cn.YIELD)
    :param bool is_test: don't plot if test
    :param dict plot_parms:
    Note: excludes alphas for which there is only one value.
    """
    if plot_parms is None:
      plot_parms = {}
    def setParms(name, value):
      """
      Updates the plot parameters.
      :parm str name: parameter name
      :parm object value:
      """
      if not name in list(plot_parms.keys()):
        plot_parms[name] = value
    #
    SIZE = 'size'
    XLABEL = "DVH Paired With Ancestral"
    YLABEL = "MMP Paired With Ancestral"
    COLOR_DICT = {"HA2": "red", "HR2": "blue", "UE3": "brown"}
    # Possible keys in plot_parms
    FIGSIZE = 'figsize'
    POINT_FONTSIZE = 'point_fontsize'
    ROTATION  = 'rotation'
    #
    setParms(FIGSIZE, (6, 4))
    setParms(POINT_FONTSIZE, 8)
    setParms(ROTATION, 5)
    #
    df = cls.makeEstimateDFS(cvsize=3, depvar=depvar)[cn.AVG]
    colors = [COLOR_DICT[Isolate.create(i).line]
              for i in df[cn.KEY_ISOLATE_DVH]]
    df.rename(columns={
        im.ALPHA_AVG: XLABEL,
        im.BETA_AVG: YLABEL,
        }, inplace=True)
    sel = [not np.isnan(x) for x in df[im.BETA_STD]]
    df = df.loc[sel].copy()
    avg = np.mean(df[cn.ESTIMATE])
    std = np.std(df[cn.ESTIMATE])
    df[SIZE] = [40*(v - avg)/std for v in df[cn.ESTIMATE]]
    df[SIZE] = df[SIZE] - np.min(df[SIZE]) + 10
    ax = df.plot(kind='scatter', x=XLABEL, y=YLABEL, c=colors, s=df[SIZE],
        title=depvar, figsize=plot_parms[FIGSIZE])
    for _, row in df.iterrows():
      label = Isolate.create(
          row[cn.KEY_ISOLATE_DVH]).getClonePairingID()
      item = ax.text(row[XLABEL]+0.12*std, row[YLABEL], label,
          rotation=plot_parms[ROTATION])
      item.set_fontsize(plot_parms[POINT_FONTSIZE])
    if not is_test:
      plt.show()


class NonParametricIsolateRegression(IsolateRegression):
  """
  Estimates values of isolates using their mean values.
  y_ijk = gamma_ij + epsilon_ijk
    where gamma_ij = mean(y_ijk)
  """

  def estimate(self, isIndex=lambda x: True, **kwargs):
    """
    Estimates alpha_ij
    :param Function isIndex:
    :return pd.DataFrame:
      cn.KEY_ISOLATE_DVH
      cn.KEY_ISOLATE_MMP
      cn.COUNT
      cn.DEPVAR
      cn.ESTIMATE
    """
    depvars = self._getDepvars(isIndex=isIndex)[cn.DEPVAR]
    estimate = np.mean(depvars)
    self.df_estimate = pd.DataFrame({
        cn.KEY_ISOLATE_DVH: [self._isolate_dvh],
        cn.KEY_ISOLATE_MMP: [self._isolate_mmp],
        cn.DEPVAR: [self._depvar],
        cn.COUNT: [self.getDataSize()],
        cn.ESTIMATE: [estimate],
        })
    return self.df_estimate


class AncestralPairingIsolateRegression(IsolateRegression):
  """
  Estimates gamma_ij as mean(y_ijk - alpha_i - beta_j) for
  i, j in the same line.
  """
  df_save = None  # Accumulated regression data

  def __init__(self, isolate_pair, depvar=cn.RATE, 
      isPermittedRow=im.ISPERMITTEDROW, 
      line=cn.LINE_HA2, leave_out=cn.KEY_CULTURE, **kwargs):
    """
    :param str depvar: Dependent variable for regression
                        cn.RATE or cn.YIELD
    :param list-of-str isolate_pair: list of isolate string
    :param function isPermittedRow: boolean function of rows
        of df_coculture
    :param str leave_out: cn.KEY_CULTURE or cn.KEY_ISOLATE
    """
    self._line = line
    self._isPermittedRow = lambda r: isPermittedRow(r)  \
        and (r[cn.LINE] == self._line)
    self._leave_out = leave_out
    super(self.__class__, self).__init__(isolate_pair,
        depvar=depvar, isPermittedRow=self._isPermittedRow)

  def _getFullIsPermittedRow(self):
    """
    Constructs a constraint
    """
    if self._leave_out == cn.KEY_ISOLATE:
      def constraint(row):
        if not self._isPermittedRow(row):
          return False
        if (row[cn.KEY_ISOLATE_DVH] == self._isolate_dvh)  \
            and (row[cn.KEY_ISOLATE_MMP] == self._isolate_mmp):
          return False
        return True
    else:
      constraint = self._isPermittedRow
    return constraint

  @classmethod
  def getKeyValues(cls, isPermittedRow=im.ISPERMITTEDROW, line=None, **kwargs):
    """
    :return list-of-list-of-str: cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP
    """
    cls._getData()
    if line is None:
      raise ValueError("Must specify line.")
    newIsPermittedRow = lambda r: isPermittedRow(r) and (r[cn.LINE] == line)
    # BUG: result returns as null. Why?
    result = ParametricIsolateRegression.getKeyValues(isPermittedRow=newIsPermittedRow,
        **kwargs)
    if len(result) == 0:
      import pdb; pdb.set_trace()
    return result

  def _saveDF(self, df):
    """
    Used as a debugging aid. Accumulates and saves dataframes.
    """
    cls = self.__class__
    df_save = df.copy()
    df_save['isolate_dvh'] = self._isolate_dvh
    df_save['isolate_mmp'] = self._isolate_mmp
    if cls.df_save is None:
      cls.df_save = df_save.copy()
    else:
      cls.df_save = pd.concat([cls.df_save, df_save])

  def estimate(self, isIndex=lambda x: True, key=None):
    """
    Estimates the mean and standard devition of values of the model constants.
    :param Function isIndex:
    :param object key: used to index models
    :return pd.DataFrame: No rows if compute estimate
      cn.KEY_ISOLATE_DVH
      cn.KEY_ISOLATE_MMP
      cn.LINE
      cn.DEPVAR
      im.ALPHA_AVG
      im.ALPHA_STD
      im.BETA_AVG
      im.BETA_STD
      COEF_DVH
      COEF_MMP
      cn.COUNT
      im.GAMMA
      cn.RSQ
      cn.ESTIMATE - estimated value for the combination of evolved isolates
    """
    def getAncestralAvgStd(key_isolate):
      results = [(r[cn.AVG], r[cn.STD])
               for _,r in df_ancestral.iterrows()
               if r[cn.KEY_ISOLATE] == key_isolate]
      if len(results) == 0:
        return None
      return results[0]
    #
    cls = self.__class__
    # Create the subsetted dataframes
    df_ancestral = cls.dfs_ancestral[self._depvar].copy()
    util.resetIndex(df_ancestral)
    df_predictor = cls.dfs_coculture[self._depvar]
    df_predictor = df_predictor[df_predictor[cn.LINE]==self._line]
    del df_predictor[cn.DEPVAR]  # Added later
    #
    def addAncestralPairingColumn(df_predictor, species):
      """
      Adds the AVG and STD values for the species.
      """
      if species == cn.SPECIES_MIX_DVH:
        avgcol = im.ALPHA_AVG
        stdcol = im.ALPHA_STD
        key_isolate = cn.KEY_ISOLATE_DVH
      else:
        avgcol = im.BETA_AVG
        stdcol = im.BETA_STD
        key_isolate = cn.KEY_ISOLATE_MMP
      sel = [Isolate.isSpecies(v, species) 
             for v in df_ancestral[cn.KEY_ISOLATE]]
      df = df_ancestral.loc[sel, :]
      df_result = df_predictor.merge(df, left_on=key_isolate,
          right_on=cn.KEY_ISOLATE, how='inner')
      del df_result[cn.KEY_ISOLATE]
      df_result.rename(columns={
          cn.AVG: avgcol,
          cn.STD: stdcol,
          }, inplace=True)
      return df_result
    #
    df_predictor = addAncestralPairingColumn(df_predictor, cn.SPECIES_MIX_DVH)
    df_predictor = addAncestralPairingColumn(df_predictor, cn.SPECIES_MIX_MMP)
    #
    def fit():
      """
      Fits a model for different predictor variables
      :return 5-tuple: coef_dvh, coef_mmp gamma, estimate, rsq
      """
      species = None
      predictors = [
          [cn.SPECIES_MIX_DVH, cn.SPECIES_MIX_MMP], 
          [cn.SPECIES_MIX_DVH],
          [cn.SPECIES_MIX_MMP],
      ]
      for predictor in predictors:
        keys = []
        if cn.SPECIES_MIX_DVH in predictor:
          keys.append(im.ALPHA_AVG)
        if cn.SPECIES_MIX_MMP in predictor:
          keys.append(im.BETA_AVG)
        if len(keys) == 0:
          return None
        X = df_predictor[keys].values
        y = df_predictor[cn.DEPVAR].tolist()
        lr = linear_model.LinearRegression(fit_intercept=True, copy_X=True)
        fit = lr.fit(X, y)
        # See if the result is co-linear
        if not util.isColinear(X):
          species = predictor
          break
      if species is None:
        raise RuntimeError("Colinear matrix.")
      #
      coef_dvh = 0.0
      coef_mmp = 0.0
      if cn.SPECIES_MIX_DVH in species:
        coef_dvh = fit.coef_[0]
      elif species == [cn.SPECIES_MIX_MMP]:
        coef_mmp = fit.coef_[0]
      if len(species) == 2:
        coef_mmp = fit.coef_[1]
      gamma = fit.intercept_
      estimates = fit.predict(X)
      rsq = fit.score(X, y)
      return coef_dvh, coef_mmp, gamma, estimates, rsq
    #
    constraint = self._getFullIsPermittedRow()
    isolate_pairs = [[r[cn.KEY_ISOLATE_DVH], r[cn.KEY_ISOLATE_MMP]]
        for _,r in cls.dfs_coculture[self._depvar].iterrows() if constraint(r)]
    isolate_pairs = util.removeDuplicatesFromList(isolate_pairs)
    # Construct the dependent variables
    dfs = []
    # Select the depvars depending on how cross validation is done
    for pair in isolate_pairs:
      constraint = lambda r: self._isPermittedRow(r) and  \
          (r[cn.KEY_ISOLATE_DVH] == pair[0]) and  \
          (r[cn.KEY_ISOLATE_MMP] == pair[1])
      if (self._isolate_dvh == pair[0]) and (self._isolate_mmp == pair[1]):
        isIndexAdj = isIndex
      else:
        isIndexAdj = lambda v: True  # Use all data if not the modelled isolate
      dfs.append(cls.getDepvarDF(self._depvar, isIndex=isIndexAdj, 
           constraint=constraint))
    df_depvar = pd.concat(dfs)
    df_predictor = df_predictor.merge(df_depvar, on=cn.KEY_CULTURE,
        how='inner')
    # Fit the model, checking for colinearity
    coef_dvh, coef_mmp, gamma, estimate, rsq = fit()
    # Construct the estimate for this isolate pair
    ancestral_dvh = getAncestralAvgStd(self._isolate_dvh)
    ancestral_mmp = getAncestralAvgStd(self._isolate_mmp)
    if (ancestral_dvh is None) or (ancestral_mmp is None):
      return pd.DataFrame()
    #
    estimate_isolate_pair = coef_dvh*ancestral_dvh[0]  \
        + coef_mmp*ancestral_mmp[0] + gamma
    row = {
        cn.LINE: [self._line],
        cn.KEY_ISOLATE_DVH: [self._isolate_dvh],
        cn.KEY_ISOLATE_MMP: [self._isolate_mmp],
        im.ALPHA_AVG: [ancestral_dvh[0]],
        im.ALPHA_STD: [ancestral_dvh[1]],
        im.BETA_AVG: [ancestral_mmp[0]],
        im.BETA_STD: [ancestral_mmp[1]],
        cn.COUNT: [self.getDataSize()],
        COEF_DVH: [coef_dvh],
        COEF_MMP: [coef_mmp],
        im.GAMMA: [gamma],
        cn.ESTIMATE: [estimate_isolate_pair],
        cn.RSQ: [rsq],
        }
    self.df_estimate = pd.DataFrame(row)
    #
    return self.df_estimate
