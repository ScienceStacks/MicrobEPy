import __init__
import helpers
import util
import constants as cn
import study_context

import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
DEPVAR = cn.RATE
MUTATION_COLUMN = cn.POSITION


class TestStudyContext(unittest.TestCase):

  def setUp(self):
    self.context = study_context.StudyContext( depvar=DEPVAR,
        mutation_column=MUTATION_COLUMN)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(self.context.depvar, DEPVAR)
    self.assertEqual(self.context.mutation_column,
        MUTATION_COLUMN)

  def testRepr(self):
    if IGNORE_TEST:
      return
    stg = str(self.context)
    for val in [DEPVAR, MUTATION_COLUMN]:
      self.assertTrue(val in stg)


class TestStudyContextIterator(unittest.TestCase):
  def setUp(self):
    self.results = []

  def testNext1(self):
    if IGNORE_TEST:
      return
    for ctx in study_context.nextStudyContext():
      self.assertTrue(ctx.depvar in cn.DEPVARS)
      self.assertTrue(ctx.mutation_column in cn.MUTATION_COLUMNS)
      self.results.append(ctx)
    self.assertEqual(len(self.results), 
        len(study_context.POSSIBLE_LINES)*  \
        len(cn.DEPVARS)*len(cn.MUTATION_COLUMNS))

  def testNext2(self):
    if IGNORE_TEST:
      return
    for ctx in study_context.nextStudyContext(
        specification=[cn.LINE, cn.DEPVAR, cn.MUTATION_COLUMN]):
      self.results.append(ctx)
    self.assertEqual(len(self.results), 
        len(cn.DEPVARS)*len(cn.MUTATION_COLUMNS)*len(
        study_context.POSSIBLE_LINES))

  def testNext3(self):
    if IGNORE_TEST:
      return
    KEY_A = "a"
    KEY_B = "b"
    KEY_C = "c"
    SPECIFICATION = {
            KEY_A: ['a:1', 'a:2'], 
            KEY_B: ['b:1', 'b:2', 'b:3'], 
            KEY_C: ['c:1'],
            }
    for ctx in  \
        study_context.nextStudyContext(specification=SPECIFICATION):
      self.results.append(ctx)
    length = 1
    for _, value in SPECIFICATION.items():
      length *= len(value)
    self.assertEqual(len(self.results), length)


if __name__ == '__main__':
    unittest.main()
