"""Abstract class for models of Isolates."""

from microbepy.common import constants as cn
from microbepy.common import util

import pandas as pd
import numpy as np

ISPERMITTEDROW = lambda r: True

ALPHA_AVG = "alpha_avg"
ALPHA_STD = "alpha_std"
BETA_AVG = "beta_avg"
BETA_STD = "beta_std"
GAMMA_AVG = "gamma_avg"
GAMMA_STD = "gamma_std"
GAMMA = "gamma"
ALPHA = "alpha"


class IsolateModel(object):
  """
  Abstract class.
  Interface used with CrossValidation. Represents a model that
  does a single estimate. The IsolateModel is responsible
  for accessing its own data.
  Key concepts:
    key_names - names of columns that constitute a key for models
    key_values - values for a key
    isIndex - a function of an int that returns a bool
              This is used to select data subsets for estimation
              and prediction
  Note that there are separate methods for estimate and predict
  since different subsets of data may be used to estimate
  parameters vs. predict values of dependent variables.
  """
  # Dictionary of dataframes keyed by RATE, YIELD
  # cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP, cn.LINE,
  # cn.KEY_CULTURE, cn.DEPVAR
  dfs_coculture = {}
  # Dictionary of dataframes for ancestral pairings. Keyed by RATE, YIELD
  # cn.KEY_ISOLATE, cn.cn.AVG, cn.cn.STD
  dfs_ancestral =  {}

  def __init__(self, key_values, isPermittedRow=ISPERMITTEDROW,
      **kwargs):
    """
    :param object model_id: identifies the model to construct
    """
    pass

  @classmethod
  def _getData(cls):
    cls._makeCocultureDFS()
    cls._makeAncestralDFS()

  @classmethod
  def resetDFS(cls):
    cls.dfs_coculture = {}
    cls.dfs_ancestral = {}

  def getDataSize(self):
    """
    Provides the number of data elements
    :return int:
    """
    raise RuntimeError("Must override!")

  @classmethod
  def makeCultureIsolateDF(cls, cocultures=None):
    """
    Creates a dataframe with the DVH and MMP isolates for the cocultures.
    :param list-of-str cocultures:
    :return pd.DataFrame: KEY_CULTURE, KEY_ISOLATE_DVH, KEY_ISOLATE_MMP
    """
    df = cls._makeCocultureDFS()[cn.RATE]
    if cocultures is None:
      cocultures = df[cn.KEY_CULTURE].tolist()
    df_result = pd.DataFrame([r for _,r in df.iterrows()
                              if r[cn.KEY_CULTURE] in cocultures])
    df_result = df_result[[cn.KEY_CULTURE, cn.KEY_ISOLATE_DVH,
        cn.KEY_ISOLATE_MMP]].copy()
    return df_result
    

  @classmethod
  def _makeAncestralDFS(cls):
    """
    Computes values in ancestral pairings.
    Creates dict-of-pd.DataFrame: keyed by cn.RATE, cn.YIELD. Columns:
       cn.KEY_ISOLATE (INDEX)
       cn.AVG
       cn.STD - standard deviation of the average
    """
    # TODO: Change when change from WT to AN
    anpat = "%s%%" % cn.LINE_AN  # Avoid python interpreting '%'
    def makeDF(depvar):
      query = '''
        select distinct key_isolate,
              key_culture as key_culture_ev, key_isolate_an, %s
          from %s,
            (select distinct key_isolate as key_isolate_an, 
             key_culture as key_culture_wt from %s) sub
         where key_culture_ev = key_culture_wt
            and key_isolate not like ("%s")
        	and key_isolate_an like("%s")
        order by key_isolate, key_isolate_an
      ''' % (depvar, cn.TABLE_CULTURE_ISOLATE_MUTATION,
          cn.TABLE_CULTURE_ISOLATE_MUTATION, anpat, anpat)
      df = util.readSQL(query)
      df = df[[cn.KEY_ISOLATE, depvar]].drop_duplicates()
      groupby = df.groupby(cn.KEY_ISOLATE)
      result = groupby.mean()
      result.rename(columns={depvar: cn.AVG}, inplace=True)
      result[cn.STD] = groupby.std()[depvar]
      result[cn.COUNT] = groupby.count()
      result[cn.STD] = [r[cn.STD]/np.sqrt(r[cn.COUNT])
                       for _,r in result.iterrows()]
      del result[cn.COUNT]
      return result
    #
    if len(cls.dfs_ancestral) == 0:
      for depvar in [cn.RATE, cn.YIELD]:
        cls.dfs_ancestral[depvar] = makeDF(depvar)
    return cls.dfs_ancestral

  @classmethod
  def _makeCocultureDFS(cls):
    """
    Computes the values of the dependent variable
    Constructs a dictionary indexed by depvar with pd.DataFrame:
      cn.KEY_CULTURE
      cn.LINE
      cn.KEY_ISOLATE_DVH
      cn.KEY_ISOLATE_MMP
      cn.DEPVAR
    """
    def makeDF(depvar):
      query = '''
      select distinct 
          key_culture, line,
          key_isolate as key_isolate_dvh, 
          key_isolate_mmp, 
          %s
        from genotype_phenotype,
            (select distinct key_isolate as key_isolate_mmp, 
                transfer as transfer_mmp,
                key_culture as key_culture_mmp from genotype_phenotype 
              where species='%s' and transfer_mmp=%d) sub
        where species_mix = '%s' 
            and transfer = %d
            and key_culture is not null 
            and key_isolate is not null
            and species = '%s' 
            and key_culture_mmp = key_culture
      ''' % (depvar, cn.SPECIES_MIX_MMP, cn.TRANSFER_1000G, 
             cn.SPECIES_MIX_BOTH, cn.TRANSFER_1000G,
             cn.SPECIES_MIX_DVH)
      df = util.readSQL(query)
      df.rename(columns={depvar: cn.DEPVAR}, inplace=True)
      return df
    #
    if len(cls.dfs_coculture) == 0:
      for depvar in [cn.RATE, cn.YIELD]:
        cls.dfs_coculture[depvar] = makeDF(depvar)
    return cls.dfs_coculture

  @classmethod
  def makeEstimateDFS(cls, cvsize=3, depvar=cn.RATE,
      isPermittedRow=ISPERMITTEDROW):
    """
    :param int cvsize: size of the cross validation set
    :param str depvar: dependent variable
    :param function isPermittedRow: True if row is permitted
        allows for filtering
    :return dict: Cross validation dataframes
    """
    cross = CrossValidation(cls, cvsize, depvar=depvar,
        isPermittedRow=isPermittedRow)
    return cross.estimate()

  @classmethod
  def getKeyValues(cls):
    """
    :return list-of-str: Names of the keys for the models
    """
    raise RuntimeError("Must override!")

  def estimate(self, isIndex=lambda x: True, **kwargs):
    """
    :param function isIndex:
       arg: int; returns bool
    :return pd.DataFrame: df_estimate
      key columns - those specified in getKeys    
      columns - names of model parameters estimated
    """
    raise RuntimeError("Must override!")

  def predict(self, isIndex=lambda x: True):
    """
    :param Function isIndex:
    :return pd.DataFrame:
      key columns - those specified in getKeys    
      cn.ESTIMATE - estimated value
      OBSERVED
    """
    raise RuntimeError("Must override!")
