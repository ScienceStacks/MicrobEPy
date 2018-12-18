import __init__
import helpers
import util
import constants as cn
from range_constraint import RangeConstraint, RangeConstraintVector

import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
LOWER = 1.0
UPPER = 2.5
IN_RANGE = 1.5
OUT_RANGE = 3.0
C1 = "c1"
C2 = "c2"
LOWER2 = 10.0
UPPER2 = 25.0
IN_RANGE2 = 15
OUT_RANGE2 = 30
CONSTRAINT_INF = RangeConstraint(lower=-np.inf, upper=np.inf)
VECTOR_INF= RangeConstraintVector({
    C1: CONSTRAINT_INF,
    C2: CONSTRAINT_INF,
    })


class TestRangeConstraint(unittest.TestCase):

  def setUp(self):
    self.constraint = RangeConstraint(lower=LOWER, upper=UPPER)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(self.constraint.lower, LOWER)

  def testRepr(self):
    if IGNORE_TEST:
      return
    stg = str(self.constraint)
    self.assertTrue("[" in stg)

  def testIsSatisfied(self):
    if IGNORE_TEST:
      return
    new_constraint = RangeConstraint(lower=np.nan, upper=np.nan)
    self.assertFalse(new_constraint.isSatisfied(OUT_RANGE))
    #
    self.assertTrue(self.constraint.isSatisfied(IN_RANGE))
    new_constraint = RangeConstraint()
    self.assertFalse(self.constraint.isSatisfied(OUT_RANGE))

  def testIsEqual(self):
    if IGNORE_TEST:
      return
    self.assertTrue(self.constraint.isEqual(self.constraint))
    constraint = RangeConstraint()
    self.assertFalse(self.constraint.isEqual(constraint))

  def testIntersection(self):
    if IGNORE_TEST:
      return
    constraint = self.constraint.intersection(self.constraint)
    self.constraint.isEqual(constraint) 

  def testIsSatisfiedRows(self):
    SIZE = 10
    LOWER = 0
    UPPER = 4.5
    df = pd.DataFrame({
        cn.VALUE: range(SIZE),
        })
    constraint = RangeConstraint(lower=LOWER, upper=UPPER)
    indices = constraint.isSatisfiedRows(df)
    df_new = df.loc[indices]
    trues = [(LOWER <= v) and (v <= UPPER) for v in df_new[cn.VALUE]]
    self.assertTrue(all(trues))


class TestRangeConstraintVector(unittest.TestCase):

  def setUp(self):
    self.vector = RangeConstraintVector({
        C1: RangeConstraint(lower=LOWER, upper=UPPER),
        C2: RangeConstraint(lower=LOWER2, upper=UPPER2),
        })

  def testConstructor(self):
    self.assertEqual(set([C1, C2]), 
        set(self.vector.constraints.keys()))

  def testRepr(self):
    stg = str(self.vector)
    for key in self.vector.constraints.keys():
      self.assertTrue(key in stg)

  def testIsSatisfiable(self):
    self.assertTrue(self.vector.isSatisfiable())
    constraints = {
        C1: RangeConstraint(lower=LOWER, upper=UPPER),
        C2: RangeConstraint(lower=np.inf, upper=-np.inf),
        }
    vector = RangeConstraintVector(constraints)
    self.assertFalse(vector.isSatisfiable())

  def testFindSatisfiedRows(self):
    df = pd.DataFrame({
        C1: [IN_RANGE, IN_RANGE, OUT_RANGE, IN_RANGE],
        C2: [IN_RANGE2, IN_RANGE2, IN_RANGE2, OUT_RANGE2],
        })
    valid_indices = set([0, 1])
    indices = self.vector.findSatisfiedRows(df)
    self.assertTrue(set(indices), valid_indices)

  def testIntersection(self):
    vector_intersection = self.vector.intersection(VECTOR_INF)
    self.assertTrue(self.vector.isEqual(vector_intersection))

  def testIsEqual(self):
    self.assertTrue(self.vector.isEqual(self.vector))
    self.assertTrue(VECTOR_INF.isEqual(VECTOR_INF))
    self.assertFalse(self.vector.isEqual(VECTOR_INF))


if __name__ == '__main__': unittest.main()
