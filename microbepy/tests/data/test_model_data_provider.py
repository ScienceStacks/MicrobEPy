from microbepy.common import constants as cn
from microbepy.common import helpers
from microbepy.common import util
from microbepy.common.isolate import Isolate
from microbepy.common.study_context import StudyContext
from microbepy.common.mutation_context import MutationContext
from microbepy.common.range_constraint  \
    import RangeConstraint, RangeConstraintVector
from microbepy.data import model_data_provider
from microbepy.data.predictor_transformer import PredictorTransformer

import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
MUTATION_CONTEXT = MutationContext(cn.RATE, cn.GGENE_ID)
    

class TestModelDataProvider(unittest.TestCase):

  def setUp(self):
    self.cls = model_data_provider.ModelDataProvider

  def testConstructor(self):
    if IGNORE_TEST:
      return
    provider = self.cls(MUTATION_CONTEXT)
    self.assertIsNotNone(self.cls.df_data)
    self.assertTrue(helpers.isValidDataFrame(self.cls.df_data,
        [cn.KEY_CULTURE, cn.LINE, cn.GENE_ID, cn.GGENE_ID,
        cn.POSITION, cn.KEY_MUTATION, cn.RATE, cn.YIELD,
        cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP],
        nan_columns=[cn.GENE_ID]))

  def testMakeXyDF(self):
    def test(constraints=None):
      provider = self.cls(MUTATION_CONTEXT, constraints=constraints)
      provider._makeXyDF()
      self.assertEqual(len(provider.df_y), len(provider.df_X))
      self.assertTrue(helpers.isValidDataFrame(
          provider.df_y,
          [cn.RATE]))
      self.assertTrue(helpers.isValidDataFrame(
          provider.df_X,
          provider.getMutations()))
      return len(provider.df_y)
    #
    len1 = test()
    constraint = lambda r: r[cn.KEY_CULTURE].count("C1") > 0
    len2 = test([constraint])
    self.assertGreater(len1, len2)
    with self.assertRaises(ValueError):
      constraint = lambda r: r[cn.KEY_CULTURE].count("C10") > 0
      test([constraint])

  def testWithRangeConstraints(self):
    if IGNORE_TEST:
      return
    def test(rc_vector=None):
      provider = self.cls(MUTATION_CONTEXT, rc_vector=rc_vector)
      provider.do()
      self.assertEqual(len(provider.df_y), len(provider.df_X))
      self.assertTrue(helpers.isValidDataFrame(
          provider.df_y,
          [cn.RATE]))
      return provider.df_y
    #
    df0 = test()
    rc_vector = RangeConstraintVector(
        {cn.RATE: RangeConstraint(lower=-1, upper=1)}
        )
    df1 = test(rc_vector)
    self.assertGreater(len(df0), len(df1))

  def makeANConstraint(self):
    return lambda r: (not Isolate.isAN(r[cn.KEY_ISOLATE_DVH]))  \
        and not (Isolate.isAN(r[cn.KEY_ISOLATE_MMP]))

  def testAggregateByIsolatePair(self):
    if IGNORE_TEST:
      return
    def test(df, tester):
      for col in df.columns:
        trues = [tester(v) for v in df[col]]
        self.assertTrue(all(trues))
    #
    provider = self.cls(MUTATION_CONTEXT)
    provider._makeXyDF()
    provider._aggregateByIsolatePair()
    test(provider.df_X,
        lambda v: np.isclose(v, 1.0) or np.isclose(v, 0.0))

  def testDo(self):
    def test(df):
      tester = lambda v: np.isclose(v, 1.0) or np.isclose(v, 0.0)
      for col in df.columns:
        trues = [tester(v) for v in df[col]]
        self.assertTrue(all(trues))
    #
    provider = self.cls(MUTATION_CONTEXT)
    provider.do()
    self.assertTrue(helpers.isValidDataFrame(provider.df_y,
        [cn.RATE]))
    self.assertTrue(helpers.isValidDataFrame(provider.df_y_std,
        [cn.RATE]))
    self.assertTrue(helpers.isValidDataFrame(provider.df_X,
        provider.df_X.columns))
    trues = [i1 == i2 for i1, i2 in 
        zip(provider.df_X.index, provider.df_y.index)]
    self.assertTrue(all(trues))
    test(provider.df_X)


  def testDoByLine(self):
    if IGNORE_TEST:
      return
    constraints = [lambda r: r[cn.LINE] == cn.LINE_UE3]
    provider1 = self.cls(MUTATION_CONTEXT, constraints=constraints)
    provider1.do()
    provider2 = self.cls(MUTATION_CONTEXT, constraints=constraints,
        is_standardize_by_line=True)
    provider2.do()
    self.assertTrue(provider1.df_y.equals(provider2.df_y))

  def testGetIsolatesFromIndices(self):
    if IGNORE_TEST:
      return
    provider = self.cls(MUTATION_CONTEXT)
    provider.do()
    isolate_dict = self.cls.getIsolatesFromIndices(provider.df_X.index)
    keys = isolate_dict.keys()
    for key in keys:
      self.assertGreater(len(isolate_dict[key]), 0)
    self.assertEqual(set(keys), 
        set([cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP]))
    self.assertEqual(len(isolate_dict[cn.KEY_ISOLATE_DVH]),
        len(isolate_dict[cn.KEY_ISOLATE_MMP]))

  def testStandardizeByLine(self):
    if IGNORE_TEST:
      return
    provider = self.cls(MUTATION_CONTEXT)
    provider.do()
    df = provider._standardizeByLine()
    self.assertTrue(helpers.isValidDataFrame(df,
        provider.df_y.columns))
    self.assertFalse(df.equals(provider.df_y))
    
    

class TestModelDataDualProvider(unittest.TestCase):

  def setUp(self):
    self.cls = model_data_provider.ModelDataDualProvider
    self.provider = self.cls(cn.GGENE_ID,
        constraints = [lambda r: r[cn.LINE] == 'HR2'])

  def testConstructor(self):
    self.assertTrue(isinstance(self.provider.df_ys, dict))

  def testDo(self):
    self.provider.do()
    for depvar in cn.DEPVARS:
      self.assertTrue(helpers.isValidDataFrame(
          self.provider.df_ys[depvar], [cn.VALUE]))
    self.assertTrue(helpers.isValidDataFrame(
          self.provider.df_X, self.provider.df_X.columns))


class TestFunctions(unittest.TestCase):

  def testMakeTransformedData(self):
    def test(provider):
      self.assertTrue(isinstance(provider, 
          model_data_provider.PredictorProvider))
      self.assertTrue(helpers.isValidDataFrame(provider.df_X,
          provider.df_X.columns.tolist()))
    #
    for data_cls in [model_data_provider.ModelDataProvider,
        model_data_provider.ModelDataDualProvider]:
      kwargs = {
          "data_cls": data_cls,
          "mutation_column": cn.GGENE_ID,
          "depvar": cn.RATE,
          }
      provider1 = model_data_provider.makeTransformedData(**kwargs)
      test(provider1)
      kwargs["constraints"] = [lambda r: r[cn.LINE] == cn.LINE_HR2]
      provider2 = model_data_provider.makeTransformedData(
          **kwargs)
      test(provider2)
      self.assertGreater(len(provider1.df_X), len(provider2.df_X))


if __name__ == '__main__':
    unittest.main()
