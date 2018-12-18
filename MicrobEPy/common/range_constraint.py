"""
Sets, tests, and manipulates interval constraints.
RangeConstraint is for a single interval.
RangeConstraintVector is for a set of intervals.
"""

import numpy as np
import constants as cn


class RangeConstraint(object):

  def __init__(self, lower=-np.inf, upper=np.inf):
    """
    :param float lower:
    :param float upper:
    """
    self.lower = lower
    self.upper = upper
    self.is_satisfiable = False

  def isSatisfiable(self):
    """
    Determines if it's possible for the constraint to be satisfied.
    :return bool:
    """
    return self.lower <= self.upper

  def isSatisfied(self, test_value):
    """
    Tests if the value lies within the range.
    :param float test_value:
    ;return bool:
    """
    result = (self.lower <= test_value)  \
        and (test_value <= self.upper)
    return result

  def isSatisfiedRows(self, df, column=cn.VALUE):
    """
    :param pd.DataFrame df:
    :param str column: Column of data to examine
    :return pd.DataFrame.Index
    """
    sel = [self.isSatisfied(v) for v in df[column]]
    return df[sel].index

  def intersection(self, range_constraint):
   """
   Calculates a RangeConstraint that is the intersection of two
   existing constraints.
   :param RangeConstraint range_constraint:
   :return RangeConstraint:
   """
   lower = max(self.lower, range_constraint.lower)
   upper = min(self.upper, range_constraint.upper)
   return RangeConstraint(lower=lower, upper=upper)

  def __repr__(self):
    return str([self.lower, self.upper])

  def isEqual(self, range_constraint):
    return (self.lower == range_constraint.lower)  \
        and (self.upper == range_constraint.upper)


class RangeConstraintVector(object):
  """Vector of range constraints."""

  def __init__(self, range_constraints):
    """
    :param dict range_constraints: value is a RangeConstraint
    """
    self.constraints = range_constraints

  def __repr__(self):
    stgs = ["%s: %s" % (str(k), str(c))
        for k, c in self.constraints.iteritems()]
    return ", ".join(stgs)

  def update(self, name, value):
    self.constraints[name] = value

  def isSatisfiable(self):
    """
    Determines if the constraints are possible to satisfy.
    :return bool:
    """
    bools = [c.isSatisfiable() for c in self.constraints.values()]
    return all(bools) 

  def findSatisfiedRows(self, df):
    """
    Locates the rows of df that are satisfied by all applicable
    RangeConstraint.
    :param pd.DataFrame df: has columns that are keys for RangeConstraints.
    :return pd.indexes: indexes of rows selected
    """
    sel = [True] * len(df)
    for key in df.keys():
      if key in self.constraints:
        values = df[key].tolist()
        sel = [s and self.constraints[key].isSatisfied(v)
            for s, v in zip(sel, values)]
    return df[sel].index

  def intersection(self, other_vector):
    """
    Creates a RangeConstraintVector that is the intersection
    of two existing RangeConstraintVectors.
    :param RangeConstraintVector other_vector:
    :return RangeConstraintVector:
    """
    keys = set(self.constraints.keys()).intersection(
        other_vector.constraints.keys())
    new_vector = {k: self.constraints[k].intersection(
        other_vector.constraints[k])
        for k in keys}
    return RangeConstraintVector(new_vector)

  def isEqual(self, other_vector):
    """
    Tests for equality between to RangeConstraintVector
    :param RangeConstraintVector other_vector:
    :return bool:
    """
    if not (self.constraints.keys() 
        == other_vector.constraints.keys()):
      return False
    tests = [self.constraints[k].isEqual(other_vector.constraints[k])
        for k in self.constraints.keys()]
    return all(tests)
