"""Tests for CorrelationStatistic"""

import unittest

import microbepy_init
import microbepy_init
import constants as cn
import helpers
import correlation_statistic as cs
import util

import copy
import math
import numpy as np
import pandas as pd


IGNORE_TEST = False
COL_A = 'a'
COL_B = 'b'
COL_C = 'c'
COL_D = 'd'
COL_E = 'e'
COL_F = 'f'


##########################################
# Classes
##########################################
class TestSetSignificanceLevel(unittest.TestCase):

  def setUp(self):
    self.df_matrix = pd.DataFrame({
        COL_A: [0, 1, 0, 1, 0, 0],
        COL_B: [0, 0, 1, 1, 1, 1],
        COL_C: [1, 1, 1, 1, 1, 1],
        COL_D: [0, 0, 0, 1, 1, 0],
        COL_E: [1, 0, 0, 0, 1, 1],
        COL_F: [0, 1, 0, 1, 0, 0],
        })
    self.ssl = cs.SetSignificanceLevel(self.df_matrix)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(self.ssl._M, len(self.df_matrix[COL_A]))

  def testCalcSet(self):
    if IGNORE_TEST:
      return
    def test(columns):
      siglvl = self.ssl.calcSet(columns)
      self.assertTrue( (siglvl <= 1.0) and (siglvl >= 0))
    #
    siglvl1 = self.ssl.calcSet([COL_B, COL_A])
    siglvl2 = self.ssl.calcSet([COL_B, COL_D])
    self.assertGreater(siglvl1, siglvl2)
    test([COL_A, COL_B])
    test([COL_A, COL_C])
    test([COL_B, COL_C])
    #
    siglvl1a = self.ssl.calcSet([COL_B, COL_A])
    # Verify that caching reproduces the original results
    self.assertTrue(np.isclose(siglvl1, siglvl1a))
    self.assertGreater(len(self.ssl._history.keys()), 0)

  def testCalcSet(self):
    if IGNORE_TEST:
      return
    def test(columns, is_sets_smaller=False):
      sets = [ [c] for c in columns]
      joint_logprob = math.log(self.ssl.calcSet(columns))
      sets_logprob = self.ssl.calcLogSets(sets)
      self.assertEqual(sets_logprob < joint_logprob, is_sets_smaller)
    #
    test([COL_A, COL_F], is_sets_smaller=False)
    test([COL_A, COL_E], is_sets_smaller=True)



##########################################
# Functions
##########################################
class TestFunctions(unittest.TestCase):

  def testNChooseR(self):
    if IGNORE_TEST:
      return
    self.assertEqual(cs.NChooseR(5, 2), 10)
    self.assertEqual(cs.NChooseR(5, 0), 1)
    self.assertEqual(cs.NChooseR(5, -1), 0)

  def testSimpleCalcProb(self):
    if IGNORE_TEST:
      return
    M = 10
    prob = cs.calcSimpleProb(M, M, M, M)
    self.assertTrue(np.isclose(prob, 1))
    prob = cs.calcSimpleProb(M, M, 1, 1)
    self.assertTrue(np.isclose(prob, 1.0))
    prob = cs.calcSimpleProb(M, 1, 1, 1)
    self.assertTrue(np.isclose(prob, 1.0/M))
    prob = cs.calcSimpleProb(M, M, 2, 1)
    self.assertTrue(np.isclose(prob, 0))

  def testCalcCumlProb(self):
    if IGNORE_TEST:
      return
    for M in [10, 20, 30]:
      for n1 in [M, M-5, M-10]:
        for n2 in [n1, int(0.8*n1), int(0.5*n1)]:
          prob = cs.calcCumlProb(M, n1, n2, 0)
          self.assertTrue(np.isclose(prob, 1.0))

  def testCalcCumlProbIntersection(self):
    if IGNORE_TEST:
      return
    def isLessEqual(v1, v2):
      if np.isclose(v1, v2):
        return True
      if v1 < v2: 
        return True
      else:
        return False
    #
    for M in [10, 20, 30]:
      for n1 in [M, int(0.8*M), int(0.5*M)]:
        for n2 in [M, int(0.8*M), int(0.5*M)]:
          prob0 = cs.calcCumlProbIntersection(M, n1, n2, 0)
          self.assertTrue(np.isclose(prob0, 1.0))
          prob1 = cs.calcCumlProbIntersection(M, n1, n2, 1)
          self.assertTrue(isLessEqual(prob1, 1.0))
          prob2 = cs.calcCumlProbIntersection(M, n1, n2, 2)
          self.assertTrue(isLessEqual(prob2, prob1))
            #print("M1: %d, M2: %d, n1: %d, n2: %d, prob(k>0): %f"
            #    %(M1, M2, n1, n2, prob))

if __name__ == '__main__':
    unittest.main()
