from microbepy.common import constants as cn
from microbepy.common import helpers
from microbepy.data.model_data_provider import ModelDataProvider
from microbepy.plot.mutation_cofraction import MutationCofraction
from microbepy.common import helpers

import numpy as np
import pandas as pd
import unittest

IGNORE_TEST = False
NUM_LINES = 9
NUM_TRANSFERS = 5
NUM_MUTATIONS = 79


class TestMutationCofraction(unittest.TestCase):

  def init(self):
    self.cofraction = MutationCofraction()

  def setUp(self):
    if IGNORE_TEST:
      return
    self.init()

  def testConstructor(self):
    if IGNORE_TEST:
      return
    o = self.cofraction
    self.assertEqual(len(o.lines), NUM_LINES)
    self.assertEqual(len(o.transfers), NUM_TRANSFERS)
    self.assertEqual(len(o.ordered_mutations), NUM_MUTATIONS)

  def testMakeLineDF(self):
    if IGNORE_TEST:
      return
    df = self.cofraction.makeLineDF(
          transfer=cn.TRANSFER_DEFAULT)
    self.assertTrue(helpers.isValidDataFrame(df, df.columns))

  def testGetFrequentMutations(self):
    if IGNORE_TEST:
      return
    mutations_2 = self.cofraction._getFrequentMutations(min_lines=2)
    mutations_3 = self.cofraction._getFrequentMutations(min_lines=3)
    self.assertGreater(len(mutations_2), len(mutations_3))

  def testMakeCoFractionDifferencedDF(self):
    if IGNORE_TEST:
      return
    df = self.cofraction.makeCofractionDifferencedDF(
        transfer=152,
        other_transfer=45)
    self.assertGreater(df.sum().sum(), 0)
    df = self.cofraction.makeCofractionDifferencedDF(
        transfer=45,
        other_transfer=45)
    self.assertEqual(df.sum().sum(), 0)

  def testMakeCoFractionDF(self):
    if IGNORE_TEST:
      return
    df = self.cofraction.makeCofractionDF(transfer=45,
        other_transfer=45, threshold_frac=0.1)
    self.assertGreater(len(df), 0)
    #
    df = self.cofraction.makeCofractionDF(transfer=45,
        other_transfer=15, threshold_frac=0.1,
        is_difference_frac=True)
    self.assertTrue(helpers.isValidDataFrame(df, df.columns))
    #
    df = self.cofraction.makeCofractionDF(transfer=15,
        other_transfer=15, threshold_frac=0.2)
    columns = df.columns
    self.assertTrue(helpers.isValidDataFrame(df, columns))
    df1 = self.cofraction.makeCofractionDF(transfer=15,
        other_transfer=15, threshold_frac=0.3)
    self.assertGreaterEqual(len(df), len(df1))


if __name__ == '__main__':
    unittest.main()
