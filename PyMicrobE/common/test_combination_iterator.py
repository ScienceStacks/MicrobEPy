import __init__
import helpers
import util
import constants as cn
from combination_iterator import CombinationIterator

import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
ELEMENTS = ['a', 'b', 'c']
MAX_SIZE = 2


class TestMutationCombination(unittest.TestCase):

  def setUp(self):
    self.combinator = CombinationIterator(ELEMENTS, MAX_SIZE)

  def testConstructor(self): 
    self.assertEqual(self.combinator._current_size, 0)

  def testNextBasic(self):
    if IGNORE_TEST:
      return
    results = []
    [results.append(c) for c in self.combinator]
    length = len(ELEMENTS)
    expected_size = length + length*(length-1)/2
    self.assertEqual(len(results), expected_size)

  def testNextExcludes(self):
    excludes = [[e] for e in ELEMENTS]
    combinator = CombinationIterator(ELEMENTS, MAX_SIZE,
        excludes=excludes)
    results = []
    [results.append(c) for c in combinator]
    length = len(ELEMENTS)
    expected_size = length*(length-1)/2
    self.assertEqual(len(results), expected_size)


if __name__ == '__main__':
    unittest.main()
