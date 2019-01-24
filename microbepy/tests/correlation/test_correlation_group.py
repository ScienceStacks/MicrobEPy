from microbepy.common import helpers
from microbepy.common import util
from microbepy.common import constants as cn
from microbepy.correlation import correlation_group as cg
from microbepy.correlation.genome_correlation import GenomeCorrelation
from microbepy.correlation import genome_correlation as gc


import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False

GENOME_CORRELATION = GenomeCorrelation()


class TestCorrelationGroup(unittest.TestCase):

  def setUp(self):
    if IGNORE_TEST:
      return
    self.gc = GENOME_CORRELATION
    self.cg = cg.CorrelationGroup(GENOME_CORRELATION,
        is_test=True,
        output_directory=cn.TEST_PATH)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(isinstance(self.cg._df_binary, pd.DataFrame))

  def testMakeKmeansGroups(self):
    if IGNORE_TEST:
      return
    def test(num_groups, is_siglvl):
      corr_group = cg.CorrelationGroup(
          GENOME_CORRELATION, is_siglvl=is_siglvl,
          is_test=True, output_directory=cn.TEST_PATH)
      df = corr_group.makeKmeansGroups(num_groups)
      self.assertEqual(set(df[cn.CATEGORICAL]),
        set(self.gc.df_binary.columns))
      valid_dict = {
          cn.GROUP: lambda x: (x >= 0) and (x < num_groups)
          }
      self.assertTrue(helpers.isValidDataFrame(df,
        [cn.CATEGORICAL, cn.GROUP],
        valid_dict=valid_dict))
    #
    test(10, False)
    test(10, True)

  def testPlotGroups(self):
    if IGNORE_TEST:
      return
    # Only a smoke test
    self.cg.plotGroups()

  def testMakeGroupedDF(self):
    if IGNORE_TEST:
      return
    def test(cluster_alg, cluster_parms):
      df = self.cg.makeGroupedDF()
      self.assertTrue(set(df.columns), 
          set(self.cg._df_corr.columns))
      self.assertGreater(len(df.index), 0)
    #
    test(cg.CATCO, cg.CATCO_DICT)
    test(cg.KMEANS, 0.3)

  def testMakeCorrelationGroupCSV(self):
    if IGNORE_TEST:
      return
    def test(cluster_alg, cluster_parms):
      df = cg.CorrelationGroup.makeCorrelationGroupCSV(
          cluster_alg=cluster_alg, cluster_parms=cluster_parms,
          output_directory=cn.TEST_PATH)
      self.assertTrue(helpers.isValidDataFrame(df,
        [cn.COUNT, cn.LINE, cn.SPECIES, cn.GGENE_ID, cn.GROUP,
        cn.GENE_DESC],
        nan_columns=[cn.GENE_DESC]))
      self.assertGreater(len(df[cn.LINE].unique()), 1)
      length_min_2 = len(df.index)
      df = cg.CorrelationGroup.makeCorrelationGroupCSV(
          output_directory=cn.TEST_PATH, min_group_size=3)
      length_min_3 = len(df.index)
      self.assertGreaterEqual(length_min_2, length_min_3)
    #
    test(cg.CATCO, {cg.MAX_SIGLVL: 0.1})
    test(cg.KMEANS, cg.FRAC_DICT)

  def testMakeCorrelationGroupDF(self):
    if IGNORE_TEST:
      return
    species = cn.SPECIES_MIX_DVH
    line_row = 'HA2'
    line_col = 'HR2'
    try:
      dfs = cg.CorrelationGroup.makeCorrelationGroupDF(species,
          line_row, line_col)
    except FileNotFoundError:
      return
    for key in list(dfs.keys()):
      df = dfs[key]
      self.assertTrue(helpers.isValidDataFrame(df, df.columns))
      for column in df.columns:
        self.assertTrue(isinstance(column, np.int64))
      for idx in df.index:
        self.assertTrue(isinstance(idx, int))

  def testMakeCorrelationGroupDF2(self):
    try:
      dfs = cg.CorrelationGroup.makeCorrelationGroupDF(cn.SPECIES_MIX_MMP, 
          'HA2', 'HR2')
    except FileNotFoundError:
      return
    self.assertIsNone(dfs)


  def testMakeCatCo(self):
    if IGNORE_TEST:
      return
    def test(species):
      constraints = [lambda r: r[cn.SPECIES] == species]
      self.cg = cg.CorrelationGroup(GENOME_CORRELATION,
          is_test=True,
          output_directory=cn.TEST_PATH)
      df_group = self.cg._makeCatCoGroups({cg.MAX_SIGLVL: 0.001})
      self.assertTrue(helpers.isValidDataFrame(df_group,
          [cn.GROUP, cn.CATEGORICAL, cn.SIGLVL]))
    #
    test(cn.SPECIES_MIX_MMP)
    test(cn.SPECIES_MIX_DVH)
    

if __name__ == '__main__':
    unittest.main()
