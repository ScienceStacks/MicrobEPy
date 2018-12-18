import __init__
import helpers
import util
import constants as cn
import mutation_context

import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
DEPVAR = cn.RATE
MUTATION_COLUMN = cn.POSITION


class TestMutationContext(unittest.TestCase):

  def setUp(self):
    self.context = mutation_context.MutationContext(
        DEPVAR, MUTATION_COLUMN)

  def testConstructor(self):
    self.assertEqual(self.context.depvar, DEPVAR)
    self.assertEqual(self.context.mutation_column,
        MUTATION_COLUMN)
    self.assertEqual(self.context.line, cn.LINE_ALL)
    self.context = mutation_context.MutationContext(
        DEPVAR, MUTATION_COLUMN, line=cn.LINE_HA2)
    self.assertEqual(self.context.line, cn.LINE_HA2)

  def testRepr(self):
    stg = str(self.context)
    for val in [DEPVAR, MUTATION_COLUMN]:
      self.assertTrue(val in stg)

  def testgetCSVName(self):
    stg = self.context.getCSVName()
    self.assertTrue("csv" in stg)


class TestMutationContextIterator(unittest.TestCase):

  def testNext(self):
    results = []
    for ctx in mutation_context.nextMutationContext():
      results.append(ctx)
    self.assertEqual(len(results), 
        len(cn.DEPVARS)*len(cn.MUTATION_COLUMNS))
    #
    results = []
    for ctx in mutation_context.nextMutationContext(
        iteration_order=[cn.LINE, cn.DEPVAR, cn.MUTATION_COLUMN]):
      results.append(ctx)
    self.assertEqual(len(results), 
        len(cn.DEPVARS)*len(cn.MUTATION_COLUMNS)*len(
        mutation_context.POSSIBLE_LINES))
    #
    results = []
    for ctx in mutation_context.nextMutationContext(
        iteration_order=[cn.LINE],
        default_values={cn.DEPVAR: cn.YIELD, 
        cn.MUTATION_COLUMN: cn.POSITION}):
      results.append(ctx)
    self.assertEqual(len(results), 
        len(mutation_context.POSSIBLE_LINES))
    trues = [(cn.YIELD in str(x)) and (cn.POSITION in str(x))
        for x in results]
    self.assertTrue(all(trues))


if __name__ == '__main__':
    unittest.main()
