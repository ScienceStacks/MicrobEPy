import microbepy_init
import constants as cn
import copy
import helpers
import util
import mutation_cooccurrence as mc
from model_data_provider import ModelDataDualProvider
from mutation_context import MutationContext
from range_constraint import RangeConstraintVector, RangeConstraint
from study_context import nextStudyContext

import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
IS_WRITE = False  # Write files with statistics


################## TEST FUNCTIONS #################
class TestCoStatistic(unittest.TestCase):

  def setUp(self):
    self.df_rate = pd.DataFrame({
      cn.ISOLATES: ['a', 'b'],
      cn.MUTATIONS: ['x', 'y'],
      cn.COUNT_MUTATIONS: [1, 1],
      cn.AVG: [1, 1],
      cn.RNG: [1, 1],
      })
    self.costatistic = mc.CoStatistic()

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(set(self.costatistic.dfs.keys()),
        set(mc.CoStatistic.SCHEMA.keys()))

  def testSet(self):
    if IGNORE_TEST:
      return
    self.costatistic.set(cn.RATE, self.df_rate)
    self.assertTrue(self.costatistic.dfs[cn.RATE].equals(self.df_rate))

  def testGet(self):
    if IGNORE_TEST:
      return
    self.costatistic.set(cn.RATE, self.df_rate)
    df_rate = self.costatistic.get(cn.RATE)
    self.assertTrue(df_rate.equals(self.df_rate))

  def testConcat(self):
    if IGNORE_TEST:
      return
    def test(df_rate):
      costatistic1 = copy.deepcopy(self.costatistic)
      costatistic1.set(cn.RATE, self.df_rate)
      costatistic2 = mc.CoStatistic()
      costatistic2.set(cn.RATE, df_rate)
      expected_length = len(costatistic1.dfs[cn.RATE])  \
          + len(df_rate)
      costatistic1.concat(costatistic2)
      self.assertEqual(len(costatistic1.dfs[cn.RATE]),
          expected_length)
      for key in [cn.MIN, cn.YIELD]:
        self.assertEqual(len(self.costatistic.dfs[key]), 0)
    #
    test(self.df_rate)
    test(pd.DataFrame())
    

################## TEST FUNCTIONS #################
class TestMutationOccurrence(unittest.TestCase):

  def doSetUp(self):
    self.provider = ModelDataDualProvider(cn.GGENE_ID,
        constraints = [lambda r: r[cn.LINE] == cn.LINE_HA2])
    self.provider.do()
    self.cooccur = mc.MutationCooccurrence(provider=self.provider)

  def setUp(self):
    if IGNORE_TEST:
      return
    self.doSetUp()

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(isinstance(self.cooccur.df_X, pd.DataFrame))

  def testFind(self):
    if IGNORE_TEST:
      return
    self.doSetUp()
    def checkEqualList(list1, list2):
      difference =  \
          set(list1).symmetric_difference(list2)
      self.assertEqual(len(difference), 0)
    #
    isolates, mutations = self.cooccur.find()
    # Check mutations
    self.assertTrue(set(mutations).issuperset(['DVU2776', 'DVU0799']))
    # Check isolates
    df_X = self.provider.df_X.copy()
    df_X = df_X.reset_index()
    expected_isolates = df_X[cn.KEY_ISOLATE_DVH].tolist()
    expected_isolates.extend(df_X[cn.KEY_ISOLATE_MMP].tolist())
    checkEqualList(expected_isolates, isolates)

  def testFindWithRangeConstraints(self):
    if IGNORE_TEST:
      return
    constraints = None
    constraints = [lambda r: r[cn.LINE] == cn.LINE_HR2]
    constraints = [lambda r: r[cn.LINE] == cn.LINE_UE3]
    cooccur = mc.MutationCooccurrence(mutation_column=cn.GGENE_ID,
        constraints=constraints)
    constraint_vector = RangeConstraintVector({
        cn.RATE: RangeConstraint(lower=-1.3, upper= -1.0),
        cn.YIELD: RangeConstraint(lower=-np.inf, upper=np.inf),
        })
    isolates, mutations = cooccur.findWithRangeConstraints(
        constraint_vector)
    self.assertTrue(len(isolates) % 2 == 0)
    self.assertGreater(len(mutations), 0)

  def testMakeCoStatistic(self):
    if IGNORE_TEST:
      return
    line = cn.LINE_HA2
    size = 5
    constraints = [lambda r: r[cn.LINE] == line]
    cooccur = mc.MutationCooccurrence(mutation_column=cn.GGENE_ID,
        constraints=constraints)
    result = cooccur.makeCoStatistic(size)
    self.assertTrue(helpers.isValidDataFrame(result.get(cn.MIN),
        mc.CoStatistic.SCHEMA[cn.MIN]))
    for depvar in cn.DEPVARS:
      self.assertTrue(helpers.isValidDataFrame(result.get(depvar),
          mc.CoStatistic.SCHEMA[depvar]))
    for df in result.values():
      trues = [v == size for v in df[cn.COUNT_ISOLATES]]
      self.assertTrue(all(trues))
    if IS_WRITE:
      for key in result.dfs.keys():
        result.get(key).to_csv("%s.csv" % key)

  def testPlot(self):
    if IGNORE_TEST:
      return
    constraints = [lambda r: r[cn.LINE] == cn.LINE_HA2]
    cooccur = mc.MutationCooccurrence(mutation_column=cn.GGENE_ID,
        constraints=constraints, is_plot=False)
    df = cooccur.plot(3)
    self.assertTrue(helpers.isValidDataFrame(df,
       [cn.RSQ, cn.SLOPE, cn.XAXIS]) )

  def init(self):
    """
    Provides common initializations.
    """
    line = cn.LINE_HA2
    size = 5
    constraints = [lambda r: r[cn.LINE] == line]
    cooccur = mc.MutationCooccurrence(mutation_column=cn.GGENE_ID,
        constraints=constraints)
    return line, size, cooccur

  def testMakeDFExplictBasic(self):
    if IGNORE_TEST:
      return
    line, size, cooccur = self.init()
    df = cooccur._makeExplicitSampleDF(size)
    expected_columns = list(cooccur.df_X.columns)
    expected_columns.extend([cn.GROUP, cn.RATE, cn.YIELD])
    self.assertTrue(helpers.isValidDataFrame(df, expected_columns))
    # Evaluate indices
    self.indexTest(df, cooccur)

  def indexTest(self, df, cooccur):
    """
    Tests that the index of a dataframe contains the correct combinations
    of DVH and MMP isolate pairs.
    :param pd.DataFrame df: index by isolate pairs
    :param mc.MutationCooccurrence cooccur:
    """
    def test(df, cooccur, idx):
      expected_size = len(df) / len(cooccur.df_X)
      values = []
      for tuple_value in df.index:
        values.append(str(tuple_value[idx]))
      ser = pd.Series(values)
      ser_count = ser.value_counts()
      trues = ser_count == expected_size
      self.assertTrue(all(trues))
    #
    test(df, cooccur, 0)
    test(df, cooccur, 1)

  def testmakeResampleSampleDF(self):
    if IGNORE_TEST:
      return
    line, size, cooccur = self.init()
    df = cooccur._makeResampleSampleDF(size, num_replications=1000)
    expected_columns = list(cooccur.df_X.columns)
    expected_columns.extend([cn.GROUP, cn.RATE, cn.YIELD])
    self.assertTrue(helpers.isValidDataFrame(df, expected_columns))

  def testMakeExplictSampleDFData(self):
    if IGNORE_TEST:
      return
    line, size, cooccur = self.init()
    explicit = cooccur.makeCoStatistic(size)
    resample = cooccur.makeCoStatistic(size, 
        num_replications=1000, is_resample=True)
    items = {"explicit": explicit, "resample": resample}
    if IS_WRITE:
      for key, costatistic in items.iteritems():
        df_rate = costatistic.get(cn.RATE)
        df_rate.to_csv("%s.csv" % key)

  def testMakeSlopeDF(self):
    if IGNORE_TEST:
      return
    df = mc.MutationCooccurrence.makeSlopeDF(
        lines=[cn.LINE_HA2], mutation_columns=[cn.GGENE_ID],
        set_sizes=[3, 4])
    self.assertTrue(helpers.isValidDataFrame(df,
        [cn.LINE, cn.RSQ, cn.SLOPE, cn.VALUE,
        cn.MUTATION_COLUMN, cn.XAXIS]))

  def testMakeLineCoStatistic(self):
    if IGNORE_TEST:
      return
    specification = {
        cn.LINE: [cn.LINE_HA2],
        cn.MUTATION_COLUMN: [cn.GGENE_ID],
        }
    rc_vector = RangeConstraintVector({
        cn.RATE: RangeConstraint(lower=-1.0, upper= 1.0),
        cn.YIELD: RangeConstraint(lower=-np.inf, upper=np.inf),
        })
    for ctx in nextStudyContext(specification):
      result = mc.MutationCooccurrence.makeLineCoStatistic(ctx,
          rc_vector=rc_vector)
      #
      for key in mc.CoStatistic.SCHEMA.keys():
        self.assertTrue(helpers.isValidDataFrame(result.get(key),
            mc.CoStatistic.SCHEMA[key]))
      #
      sizes = [len(df) for df in result.values()]
      self.assertEqual(len(set(sizes)), 1)
            

if __name__ == '__main__':
    unittest.main()
