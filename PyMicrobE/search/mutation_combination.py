"""Finds significant mutation combinations."""

import constants as cn
from group_significance_level import GroupSignificanceLevel
from model_data_provider import ModelDataProvider
from combination_iterator import CombinationIterator
import util
import util_data as ud

import itertools
import numpy as np
import pandas as pd



##############################################
# Class
##############################################
class MutationCombination(object):

  def __init__(self, mutation_context,
      transform_type=cn.TRANSFORM_LOW_FREQUENCY_ISOLATES,
      constraints=None, lines=None):
    """
    :param MutationContext mutation_context:
    :param str transform_type:
    :param list-BooleanFunction constraints: no line constraints
    :param list-str lines: Lines for which analysis is done
    """
    self._lines = util.setNoneList(lines)
    self._context = mutation_context
    self._transform_type = transform_type
    self._constraints = util.setNoneList(constraints)

  @classmethod
  def isMock(cls):
    """
    Used for dependency injection.
    """
    return False

  @classmethod
  def getFilePrefix(cls):
    """
    Used for dependency injection.
    """
    return ""

  def do(self, max_combination, is_tstat=True,
      is_resample=True,
      excludes=None,
      num_combinations=None,
      lines=None):
    """
    Searches combinations of mutations and reports their statistical
    significance.
    :param int max_combination: max mutations in a combination
    :param bool is_tstat: report t statistic
    :param bool is_resample: report resample statistic
    :param list-object excludes: combinations to exclude
    :param list-str lines: lines for which analysis is done
    :param int num_combinations: maximum combinations computed
                                  all if None
    :return pd.DataFrame, pd.DataFrame: DFs for extremas df_min, df_max 
        cn.MUTATIONS, 
        cn.SL_TSTAT
        cn.SL_RESAMPLE
        cn.VALUE - value for the extrema
        cn.COUNT - number of values in the extrema
        cn.GROUP - group for the extrema in binary
        cn.LINE
    """
    def assignDF(sl, mutation_combination, line):
      df = pd.DataFrame({
          cn.AVG: [sl.avg],
          cn.GROUP: [sl.group],
          cn.COUNT: [sl.count],
          cn.MUTATIONS: [mutation_combination],
          cn.LINE: [line],
          cn.SL_TSTAT: [sl.sl_tstat],
          cn.SL_RESAMPLE: [sl.sl_resample],
          })
      return df
    #
    excludes = util.setNoneList(excludes)
    # Get the mutations to consider
    df_min = pd.DataFrame()
    df_max = pd.DataFrame()
    lines = util.setNoneList(lines)
    if len(lines) == 0:
      lines = [cn.LINE_ALL]
    #
    combination_count = 0
    done = False
    for line in lines:
      if done:
        break
      constraints = list(self._constraints)
      if line != cn.LINE_ALL:
        constraints.append(lambda r: r[cn.LINE] == line)
      provider = ModelDataProvider(self._context, constraints=constraints)
      provider.do(transform_type=self._transform_type)
      mutations = provider.getMutations()
      combinator = CombinationIterator(mutations, max_combination,
          excludes=excludes)
      count = 0
      for vals in combinator:
        if num_combinations is not None:
          if count > num_combinations:
            done = True
            break
        count += 1
        mutation_combination = [str(m) for m in vals]
        group_sl = GroupSignificanceLevel(provider, 
            mutation_combination)
        if is_tstat:
          group_sl.calcMinMaxSLTStat()
        if is_resample:
          group_sl.calcMinMaxSLResample()
        sl_min = group_sl.sl_min
        sl_max = group_sl.sl_max
        avg = sl_min.avg
        group = sl_min.group
        count = sl_min.count
        df = assignDF(sl_min, mutation_combination, line)
        df_min = df_min.append(df, ignore_index=True)
        df = assignDF(sl_max, mutation_combination, line)
        df_max = df_max.append(df, ignore_index=True)
    return df_min, df_max
