"""
Specifies the context for a mutation or mutation combination.
A context includes:
  dependent variable
  mutation granularity
  extrema (optional), either min or max
"""

from microbepy.common import constants as cn
from microbepy.common import util

import itertools

POSSIBLE_LINES = [cn.LINE_HA2, cn.LINE_HR2, cn.LINE_UE3]


class MutationContext(object):

  def __init__(self, depvar, mutation_column, line=cn.LINE_ALL):
    if not depvar in cn.DEPVARS:
      raise ValueError("Invalid dependent: %s" % depvar)
    self.depvar = depvar
    if not mutation_column in cn.MUTATION_COLUMNS:
      raise ValueError("Invalid mutation_column: %s" 
          % mutation_column)
    self.mutation_column = mutation_column
    if not ((line in POSSIBLE_LINES) or (line == cn.LINE_ALL)):
      raise ValueError("Invalid line: %s" % line)
    self.line= line

  def __repr__(self):
    return  "%s-%s-%s" % (
        self.depvar, str(self.mutation_column), self.line)

  def getCSVName(self):
    """
    :return str: Name of the CSV file for this context
    """
    return "%s.csv" % self.__repr__()


###############################################################
def nextMutationContext(
    iteration_order=[cn.DEPVAR, cn.MUTATION_COLUMN],
    default_values={cn.DEPVAR: cn.RATE, 
    cn.MUTATION_COLUMN: cn.GGENE_ID, cn.LINE: cn.LINE_ALL}):
  """
  Generator returns the next MutationContext
  :param list-str iteration_order: Order in which values are returned
  :return MutationContext:
  """
  mega_list = []
  index_dict = {}
  pos = 0
  # Set up the order in which values will change
  for item in iteration_order:
    if item == cn.DEPVAR:
      mega_list.append(cn.DEPVARS)
      index_dict[cn.DEPVAR] = pos
    elif item == cn.MUTATION_COLUMN:
      mega_list.append(cn.MUTATION_COLUMNS)
      index_dict[cn.MUTATION_COLUMN] = pos
    elif item == cn.LINE:
      mega_list.append(POSSIBLE_LINES)
      index_dict[cn.LINE] = pos
    else:
      raise ValueError("Invalid iteration type %s" % item)
    pos += 1
  # Iterate across values
  def getValue(name):
    if name in index_dict.keys():
      value = combination[index_dict[name]]
    else:
      value = default_values[name]
    return value
  #
  iterator = itertools.product(*mega_list)
  for combination in iterator:
    depvar = getValue(cn.DEPVAR)
    line = getValue(cn.LINE)
    mutation_column = getValue(cn.MUTATION_COLUMN) 
    yield MutationContext(depvar, mutation_column, line=line)
