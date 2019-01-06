"""
Specifies the context for a study.
A study consists of values of one or more variables that are
changed such as the dependent variable, mutation granularity,
and line.
"""

import constants as cn
import util

import itertools

# Default names
# cn.DEPVAR - dependent variable
# cn.LINE - line
# cn.MUATION_COLUMN

# Default values
POSSIBLE_LINES = [cn.LINE_HA2, cn.LINE_HR2, cn.LINE_UE3]
DEFAULT_SPECIFICATION = {
    cn.DEPVAR: cn.DEPVARS,
    cn.MUTATION_COLUMN: cn.MUTATION_COLUMNS, 
    cn.LINE: POSSIBLE_LINES,
    }


class StudyContext(dict):
  """Extends dictionary with a values based representation."""

  def __init__(self, **kwargs):
    """
    :param dict kwargs: key is variable, value is its value
    """
    self.kwargs = kwargs
    # Create instance variables
    for key, value in kwargs.items():
      setattr(self, key, value)

  def __repr__(self):
    values = [str(v) for v in list(self.kwargs.values())]
    return "--".join(values)


###############################################################
def nextStudyContext(specification=DEFAULT_SPECIFICATION):
  """
  Generator returns the next StudyContext
  :param dict or list specification: 
      dict:
        key is variable name
        value is list
      list: sequence of known keys with default values
  :return StudyContext:
  """
  if isinstance(specification, list):
    new_specification = {}
    for key in specification:
      if key in list(DEFAULT_SPECIFICATION.keys()):
        new_specification[key] = DEFAULT_SPECIFICATION[key]
      else:
        raise ValueError("Invalid specification. No defaults for %s: " % key)
  else:
    new_specification = specification
  specification = new_specification
  #
  iterator = itertools.product(*list(specification.values()))
  for combination in iterator:
    kwargs = {k: v for k,v in 
        zip(list(specification.keys()), combination)}
    yield StudyContext(**kwargs)
