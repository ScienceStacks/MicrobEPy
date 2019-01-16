from microbepy.common import constants as cn
from microbepy.common import helpers
from microbepy.common.mutation_context import MutationContext
from microbepy.common import util
from microbepy.data import util_data as ud
from microbepy.data.model_data_provider import ModelDataProvider

import numpy as np
import pandas as pd
import scipy.stats as stats
import unittest


IGNORE_TEST = True

MUTATION_CONTEXT = MutationContext(cn.RATE, cn.GENE_ID)

COL_A = 'a'
COL_B = 'b'
COL_Y = 'y'
SIZE = 8
DF_X = pd.DataFrame({
    COL_A: [0, 0, 1, 1],
    COL_B: [1, 1, 0, 0],
    cn.KEY_ISOLATE_DVH: 
        ['HA2.152.01.01.D.CI', 'HA2.152.01.02.D.CI',
        'HA2.152.01.03.D.CI', 'HA2.152.01.04.D.CI'],
    cn.KEY_ISOLATE_MMP: 
        ['HA2.152.01.01.M.CI', 'HA2.152.01.02.M.CI',
        'HA2.152.01.03.M.CI', 'HA2.152.01.04.M.CI'],
    })
DF_X = DF_X.set_index([cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP])
DF_Y = pd.DataFrame({
    COL_Y: [1, 3, 5, 7],
    })

################## TEST FUNCTIONS #################
class TestFunctions(unittest.TestCase):

  def testMakeBaseDF(self):
    if IGNORE_TEST:
      return
    df_data = ud.makeCultureIsolateMutationDF()
    self.assertTrue(helpers.isValidDataFrame(df_data,
      [cn.LINE, cn.KEY_CULTURE, cn.KEY_MUTATION, cn.GGENE_ID,
      cn.KEY_ISOLATE_MMP, cn.KEY_ISOLATE_DVH,
      cn.GENE_ID, cn.RATE, cn.YIELD],
      nan_columns=[cn.GENE_ID]))

  def testMakeGausianMixture(self):
    if IGNORE_TEST:
      return
    num_rows = 10
    num_extracols = 2
    means = [1, 10]
    stds = [1, 1]
    df_X, df_y = ud.makeGausianMixture(means, stds, num_rows, 
        num_extracols)
    self.assertEqual(len(df_X), len(df_y))
    self.assertEqual(len(df_X), num_rows)
    predictor_columns = df_X.columns[0:len(means)]
    trues = ["X_" in c for c in predictor_columns]
    self.assertTrue(all(trues))
    self.assertEqual(len(df_X.columns), len(means) + num_extracols)

  def testMakeGroups(self):
    if IGNORE_TEST:
      return
    SIZE = 10
    columns = ['a', 'b', 'c']
    df = ud.makeGroups(SIZE, columns)
    self.assertEqual(len(df), SIZE)
    self.assertEqual(set(df.columns), set(columns))

  def testMakeNoisyOr(self):
    if IGNORE_TEST:
      return
    def test(num_rows, num_predictors, noise_prob):
      predictor_probs = np.repeat(
           1.0/num_predictors, num_predictors)
      df_X, df_y,score = ud.makeNoisyOr(
          predictor_probs, num_rows, noise_prob)
      self.assertEqual(len(df_X), len(df_y))
      self.assertEqual(len(df_X), num_rows)
      predictor_columns = df_X.columns[0:len(predictor_probs)]
      trues = ["P_" in c for c in predictor_columns]
      self.assertTrue(all(trues))
      self.assertEqual(len(predictor_columns), num_predictors)
    #
    test(1000, 5, 0.2)
    test(10, 2, 0.2)
    test(10, 2, 1.0)
    test(10, 2, 0)

  def testMakeNoisyBin(self):
    if IGNORE_TEST:
      return
    def test(num_rows, noise_prob):
      df_X, df_y = ud.makeNoisyBin(num_rows, noise_prob)
      predictor_columns = df_X.columns
      self.assertGreaterEqual(len(predictor_columns), np.log2(num_rows))
      trues = ["P_" in c for c in predictor_columns]
      self.assertTrue(all(trues))
      if num_rows > 100:
        # Must have enough samples to do these tests
        num_ones = df_y.sum()[cn.VALUE]
        self.assertLess(num_ones, len(df_y)*(0.75))
        self.assertGreater(num_ones, len(df_y)*(0.3))
        expected = num_rows*(1 + noise_prob)
        self.assertLess(len(df_X), expected*(1.05) )
    #
    test(20, 0.2)
    test(200, 0.2)
    test(20, 0)
    test(200, 0)

  def testMakeClassificationData(self):
    if IGNORE_TEST:
      return
    SIZE =100
    DATA = range(SIZE)
    df_X = pd.DataFrame({
      'a': DATA,
      'b': DATA,
      })
    df_y = pd.DataFrame({cn.VALUE: DATA})
    def test(percentile):
      new_df_X = df_X.copy()
      new_df_y = df_y.copy()
      ud.makeClassificationData(new_df_X, new_df_y, cn.VALUE,
          percentile)
      expected_count = int(2*percentile*len(df_y))
      self.assertLess(len(new_df_y), expected_count + 1)
    #
    test(50)
    test(10)

  def testMakeNoisyBin2(self):
    if IGNORE_TEST:
      return
    num_rows = 1000
    ONE_PROBS = [0.3, 0.5, 0.8]
    NOISE_PROBS = [0.1, 0.5, 0.7]
    NOISE_PROBS = [0.5, 0.7]
    for one_prob in ONE_PROBS:
      # Probability of a 1
      predicate=lambda v: stats.bernoulli.rvs(one_prob, size=1)[0]
      for noise_prob in NOISE_PROBS:
        _, df_y = ud.makeNoisyBin(num_rows, noise_prob,
          predicate=predicate)
        fract_one = (1.0*df_y[cn.VALUE].sum())/len(df_y)
        self.assertLess(abs(fract_one - one_prob), 0.1)

  def testMakeMutationGroupDF(self):
    if IGNORE_TEST:
      return
    df = ud.makeMutationGroupDF(DF_X, DF_Y, [COL_A, COL_B])
    self.assertTrue(helpers.isValidDataFrame(df,
        [cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP,
        cn.GROUP, cn.LINE, cn.DEPVAR]))

  def testMakeMutationGroupDFWithNumbers(self):
    if IGNORE_TEST:
       return
    COL_1 = 1
    COL_2 = 2
    df_X = DF_X.rename(columns={
        COL_A: COL_1,
        COL_B: COL_2,
        })
    df = ud.makeMutationGroupDF(df_X, DF_Y, [COL_1, COL_2])
    self.assertTrue(helpers.isValidDataFrame(df,
        [cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP,
        cn.GROUP, cn.LINE, cn.DEPVAR]))

  def testGeneratePhenotypeData(self):
    if IGNORE_TEST:
       return
    NUM_GROUP = 4
    VALUES_PER_GROUP = 5
    DATA_SIZE = NUM_GROUP*VALUES_PER_GROUP
    DATA = range(DATA_SIZE)
    GROUPS = np.repeat(VALUES_PER_GROUP, NUM_GROUP).tolist()
    NUM_REPLICATION = 50
    df = ud.generatePhenotypeData(DATA,
        GROUPS, NUM_REPLICATION)
    # Verify basics
    self.assertEqual(NUM_REPLICATION*DATA_SIZE, len(df))
    self.assertTrue(helpers.isValidDataFrame(df,
        [cn.REPLICATION, cn.GROUP, cn.VALUE]))
    # Verify that group means differe
    dfg = df.groupby([cn.GROUP]).mean()
    self.assertGreater(np.std(dfg[cn.VALUE]), 0)
    # Verify that same data are used in all replications
    dfg = df.groupby([cn.REPLICATION]).mean()
    self.assertEqual(np.std(dfg[cn.VALUE]), 0)

  def testGeneratePhenotypeDataForMutations(self):
    if IGNORE_TEST:
       return
    MUTATIONS = ['DVU1092', 'MMP1612']
    provider = ModelDataProvider(MUTATION_CONTEXT)
    provider.do()
    df = ud.generatePhenotypeDataForMutations(provider.df_X,
        provider.df_y, MUTATIONS, 1)
    self.assertTrue(helpers.isValidDataFrame(df,
        [cn.REPLICATION, cn.GROUP, cn.VALUE]))
    self.assertEqual(df[cn.GROUP].unique().tolist()[0], 0)
    self.assertEqual(len(df[cn.REPLICATION].unique()), 1)

  def testStripMutationGroup(self):
    if IGNORE_TEST:
       return
    def test(group, expected):
      result = ud.stripMutationGroup(group)
      trues = [r == e for r, e in zip(result, expected)]
      self.assertTrue(all(trues))
    #
    result = test(["aaa"], ["aaa"])
    result = test(["aaa" + ud.MUTATION_GROUP_STRING + "ccc"], ["aaa"])
    result = test(["aaa", "bbb"], ["aaa", "bbb"])
    result = test(["aaa" + ud.MUTATION_GROUP_STRING + "ccc", "bbb"], 
        ["aaa", "bbb"])

  def testMakeIsolateData(self):
    if IGNORE_TEST:
       return
    df = ud.makeIsolateData()
    self.assertTrue(helpers.isValidDataFrame(df, df.columns,
        nan_columns=df.columns))

  def testMakeStandardizedResidualsForPairings(self):
    if IGNORE_TEST:
       return
    df = ud.makeStandardizedResidualsForPairings()
    self.assertTrue(helpers.isValidDataFrame(df,
        [cn.KEY_CULTURE, ud.RATE_RES, ud.YIELD_RES]))

  def testFilterOutlierCultures(self):
    df_initial = ud.makeCultureIsolateMutationDF()
    df = ud.filterOutlierCultures(df_initial)
    self.assertGreater(len(df_initial), len(df))
    self.assertTrue(helpers.isValidDataFrame(df, df_initial.columns,
        nan_columns=df_initial.columns))

if __name__ == '__main__':
    unittest.main()
