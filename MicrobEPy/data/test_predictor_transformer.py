import __init__
import helpers
import util
import genome_model as gm
from model_data_provider import ModelDataProvider
from predictor_transformer import PredictorTransformer
from mutation_context import MutationContext

import constants as cn
import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
    

class TestPredictorTransformer(unittest.TestCase):

  def setUp(self):
    mutation_col = cn.GGENE_ID
    provider = ModelDataProvider(
        MutationContext(cn.RATE, mutation_col))
    provider.do()
    self.df_X = provider.df_X
    self.transformer = PredictorTransformer(self.df_X, mutation_col) 

  def testMakeCultureIsolate(self):
    if IGNORE_TEST:
      return
    self.assertTrue(helpers.isValidDataFrame(
        PredictorTransformer.df_culture_isolate,
        [cn.KEY_CULTURE, cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_DVH]))

  def testGetLowOccurrenceMutations(self):
    if IGNORE_TEST:
      return
    mutations = self.transformer._getLowLineOccurrenceMutations(10)
    self.assertTrue(isinstance(mutations, list))
    self.assertGreater(len(mutations), 0)

  def testTransformDefault(self):
    if IGNORE_TEST:
      return
    self.transformer.transformDefault()
    self.assertGreater(len(self.df_X.columns),
        len(self.transformer.df_X.columns))
    found = False

  def testMakeMutationGroups(self):
    if IGNORE_TEST:
      return
    self.transformer.makeMutationGroups()
    found = False
    for col in self.transformer.df_X.columns:
      if '--' in col:  # Merged columns
        found = True
        break
    self.assertTrue(found)

  def testFilterForMultipleLines(self):
    if IGNORE_TEST:
      return
    self.transformer.makeMutationGroups()
    count1 = len(self.transformer.df_X.columns)
    self.transformer.filterForMultipleLines(10)
    count2 = len(self.transformer.df_X.columns)
    self.assertGreater(count1, count2)

  def testFilterHighFrequencyMutations(self):
    self.transformer.makeMutationGroups()
    count1 = len(self.transformer.df_X.columns)
    self.transformer.filterHighFrequencyMutations()
    count2 = len(self.transformer.df_X.columns)
    self.assertGreater(count1, count2)

  def testFilterLowFrequencyColumns(self):
    if IGNORE_TEST:
      return
    count1 = len(self.transformer.df_X.columns)
    self.transformer.filterLowFrequencyColumns(300,
        min_shared_isolates=0)
    count2 = len(self.transformer.df_X.columns)
    self.assertEqual(count1, count2)
    #
    self.transformer.filterLowFrequencyColumns(300,
        min_shared_isolates=2)
    count2a = len(self.transformer.df_X.columns)
    self.assertLess(count2a, count2)
    #
    max_columns = 5
    self.transformer = PredictorTransformer(self.df_X, 
        cn.GGENE_ID)
    self.transformer.filterLowFrequencyColumns(max_columns)
    count3 = len(self.transformer.df_X.columns)
    self.assertGreaterEqual(max_columns, count3)

  def testTransformLowFrequencyColumns(self):
    if IGNORE_TEST:
      return
    count1 = len(self.transformer.df_X.columns)
    self.transformer.transformLowFrequencyColumns(3)
    count3 = len(self.transformer.df_X.columns)
    self.assertGreater(count1, count3)

if __name__ == '__main__':
    unittest.main()
