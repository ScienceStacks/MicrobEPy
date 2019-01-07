"""
Provide Statistics and Tests for the Significance of Correlations

Consider a population of size M where each member of the population
(referred to as an item) has the same  binary variables inferred by 
the presence (true) or absence (false) of a categorical value for
that item.  We are interested in the relationships between
binary variables. Specifically given a pair of binary variables,
does each item in the sample satisfy an exclusive-or 
(either both true or both false) for the binary variables?

Consider a pair of binary variables. The first is true for
m1 items (m1 <= M) and the second for m2 items. Further,
they are both true in k <= min(m1, m2).

Our null hypothesis is that the probability that a binary variable
is true is independent and identically distributed for the
items. Under our hypothesis, binary variable 1 is true with
probability m1/M for an item. Similarly, variable 2 is true
with probability m2/M. So, the probbility of k co-occurrences of
true of the two binary variables in a binomial in k and
p = m1/M*m2/M.
"""

import microbepy_init
import constants as cn

import math
import operator as op
import pandas as pd
import numpy as np
from scipy.stats import binom

LOG_OF_SMALL_VALUE = -100


#####################################
# Classes
#####################################
class SetSignificanceLevel(object):
  """
  Computes the significance level of the co-occurrence of attributes
  using a binomial distribution as the null hypothesis.
  This class assumes that the provided dataframe is not modified
  during the lifetime of an instantiated SetSignificanceLevel object.
  """

  def __init__(self, df_matrix):
    """
    :param pd.DataFrame df: rows are attributes, columns are instances,
        values are bools indicating the presence or absence of an 
        attribute.
    """
    self._df_matrix = df_matrix
    self._M = len(self._df_matrix.index)
    self._history = {}  # Diciontary of calculations

  def _calcIndProb(self, columns):
    """
    Computes the probability of co-occurrences of the binary variable
    if they are independent.
    :parm list-of-str columns:
    :return float:
    """
    df = self._df_matrix[columns]
    series = 1.0*df.sum(axis=0) / len(df.index)
    return series.product()

  def calcSet(self, columns):
    """
    Computes the significance level of the co-occurrence of
    a set of attributes (columns)
    under the null distribution that attributes are independent.
    This is a 1-sided test; that is, it does not test for very
    low values of co-occurrences.
    :param list-of-str columns: Columns that constitute the set
    :return float: significance level of co-occurrence (the
        tail of a binomial distribution)
    """
    column_list = list(columns)
    column_list.sort()
    column_hash = '--'.join(column_list)
    if column_hash in list(self._history.keys()):
      return self._history[column_hash]
    df_matrix = pd.DataFrame(self._df_matrix[column_list])
    df_transpose = df_matrix.transpose()
    k_count = sum([df_transpose[s].product() 
             for s in df_transpose.columns])
    # Obtain the tail of the distribution including the probability
    # of the number of counts observed
    if k_count > 0:
      joint_prob = self._calcIndProb(column_list)
      result = 1 - binom.cdf(k_count-1, self._M, joint_prob)
    else:
      result = 1.0
    result = max(min(result, 1.0), 0.0)  # Handle floating point
    self._history[column_hash] = result
    return result

  def calcLogSets(self, sets):
    """
    Calculates the log of the product of the significance levels.
    :param list-of-set-of-columns:
    :return float:
    """
    try:
      result = sum([math.log(self.calcSet(s)) for s in sets])
    except ValeError:
      # Encountered a 0
      result = LOG_OF_SMALL_VALUE
    return result

         

#####################################
# Functions
#####################################
def multlist(values):
  result = 1
  for value in values:
    result *= value
  return result

def NChooseR(n, r):
  if r < 0:
    return 0
  if r == 0:
    return 1
  if n < 1:
    return 0
  r = min(r, n-r)
  if r == 0: return 1
  numer = multlist(range(n, n-r, -1))
  denom = multlist(range(1, r+1))
  return numer//denom

def calcSimpleProb(M, n_1, n_2, k):
  """
  Calculates the probability of the co-occurrence of two binary
  variables in M samples
  :param int M:
  :param int n_1:
  :param int n_2:
  :param int k:
  :return float:
  """
  m1 = max(n_1, n_2)
  m2 = min(n_1, n_2)
  #
  def calc(k):
    return NChooseR(M, k) * NChooseR(M - k, m1 - k) * NChooseR(M - m1, m2 - k)
  #
  if m1 > M or m2 > M:
    return 0
  if k > m1 or k > m2:
    return 0
  if m1 + m2 > M:
    if k < m1 + m2 - M:
      return 0
  if k > min(m1, m2):
    return 0
  #
  denominator = NChooseR(M, m1) * NChooseR(M, m2)
  numerator = float(calc(k))
  result = numerator / denominator
  return result

def calcCumlProb(M, m1, m2, k):
  m_min = min(m1, m2)
  result = sum([calcSimpleProb(M, m1, m2, kk) for
                kk in range(k, m_min + 1)])
  return result

def calcCumlProbIntersection(M, m1, m2, k):
  """
  Calculates the cumulative probability of at least k itmes
  (out of a total of M) have a true value for two randomly
  assigned binary random variables, where the first is true
  in m1 items and the second is true in m2 items.
  This is an approximation that is an upper bound.
  :param int M: Number of items
  :param int m1: Number of items for which the first binary variable
      is true
  :param int m2: Number of items for which the second binary variable
      is true
  :param int k:  Number of items in which both are true
  :return float: probability
  """
  if (m1 > M) or (m2 > M):
    import pdb; pdb.set_trace()
    raise ValueError("Invalid parameters")
  # Make nmin the min and nmax the max
  if m1 < m2:
    m_min = m1
    m_max = m2
  else:
    m_min = m2
    m_max = m1
  #
  if k > m_min:
    return 0.0
  if k == 0:
    return 1.0
  # Binomial probability
  prob = float(m_min)/M * float(m_max)/M
  result = 1 - binom.cdf(k - 1, m_min, prob)
  return result
