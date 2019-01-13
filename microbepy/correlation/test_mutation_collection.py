import microbepy
from microbepy.common import constants as cn
from microbepy.correlation import group_collection
from microbepy.common import helpers
from microbepy.correlation.mutation_collection  \
    import MutationCollection
from microbepy.data import model_data_provider
from microbepy.common import util

import copy
import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
IS_WRITE = False  # Write files with statistics
LINE = "HA2"
LINE = "UE3"
IS_PLOT = False


class TestMutationCollection(unittest.TestCase):

  def setUp(self):
    constraints = [lambda r: r[cn.LINE] == LINE]
    self.provider = model_data_provider.makeTransformedData(
        mutation_column=cn.GGENE_ID,
        constraints=constraints)
    groups = MutationCollection.makeGroups(
        self.provider)
    self.collection = MutationCollection(groups, is_plot=IS_PLOT)
    self.collection.setGroupLabels(prefix=LINE)

  def testMakeGroups(self):
    if IGNORE_TEST:
      return
    collection = MutationCollection.makeGroups(
        self.provider)
    self.assertGreater(len(collection), 0)
    self.assertTrue(isinstance(collection[0], 
        microbepy.common.group_collection.Group))

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertGreater(len(self.collection.groups), 0)

  def testCalcSL(self):
    if IGNORE_TEST:
      return
    prob = MutationCollection.calcSL(10, 2, 10)
    self.assertTrue(np.isclose(prob, 1.0))
    prob1 = MutationCollection.calcSL(10, 2, 5)
    prob2 = MutationCollection.calcSL(10, 4, 5)
    self.assertGreater(prob1, prob2)

  def testMakeMutationCollectionForLine(self):
    if IGNORE_TEST:
      return
    collection = MutationCollection.makeMutationCollectionForLine(
        line=LINE, is_plot=IS_PLOT)
    self.assertTrue(collection.equals(self.collection))
    collection_all = MutationCollection.makeMutationCollectionForLine(
        is_plot=IS_PLOT)
    self.assertGreater(len(collection_all.groups), 
        len(collection.groups))
    
  def testPlot(self):
    if IGNORE_TEST:
      return
    # Smoke tests
    self.collection.plot()
    collection = MutationCollection.makeMutationCollectionForLine(
        line=cn.LINE_HR2, is_plot=IS_PLOT)
    collection.plot()
    collection = MutationCollection.makeMutationCollectionForLine(
        is_plot=IS_PLOT)
    collection.plot()
    collection = MutationCollection.makeMutationCollectionForLine(
        species=cn.SPECIES_MIX_MMP, is_plot=IS_PLOT)
    collection.plot()
    

if __name__ == '__main__':
    unittest.main()
