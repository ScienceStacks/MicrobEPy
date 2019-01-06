"""Tests for GenomeCorrelation."""

import unittest

import __init__
import constants as cn
import genome_correlation
from genome_correlation import GenomeCorrelation
import genome_correlation as gc
import helpers
import numpy as np
import pandas as pd
import util


IGNORE_TEST = False
DF15 = pd.DataFrame({
    'HA2': [1, 1],
    'HA3': [0, 0],
    'HE2': [0, 0],
    'HE3': [1, 0],
    'HR1': [0, 0],
    'HR2': [0, 0],
    'HS3': [1, 0],
    'UA2': [0, 0],
    'UA3': [0, 1],
    'UE2': [0, 0],
    'UE3': [1, 0],
    'UR1': [0, 0],
    'US1': [0, 0],
    })
DF15.index = ['MMP0234', 'DVU0799']
DF15 = DF15.transpose()
DF152 = pd.DataFrame({
    'HA2': [1, 1],
    'HA3': [0, 1],
    'HE2': [0, 1],
    'HE3': [0, 1],
    'HR1': [0, 1],
    'HR2': [0, 1],
    'HS3': [0, 0],
    'UA2': [0, 1],
    'UA3': [1, 1],
    'UE2': [1, 1],
    'UE3': [0, 1],
    'UR1': [0, 1],
    'US1': [0, 1],
    })
DF152.index = DF15.columns.tolist()
DF152 = DF152.transpose()


class TestGenomeCorrelation(unittest.TestCase):

  def setUp(self):
    self.corr = GenomeCorrelation(
        instance_name=cn.KEY_ISOLATE,
        categorical_name=cn.GENE_ID,
        is_test=True)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(helpers.isValidDataFrame(
        self.corr.df_binary, self.corr.df_binary.columns))

  def testMakeBaseDF(self):
    if IGNORE_TEST:
      return
    columns = [cn.KEY_MUTATION, cn.KEY_ISOLATE, cn.COUNT,
        cn.GENE_ID, cn.COMMUNITY, cn.SPECIES, cn.LINE]
    self.assertTrue(helpers.isValidDataFrame(
        self.corr.df_base, columns,
        nan_columns=[cn.GENE_ID]))

  def testMakeCorrelationDF(self):
    if IGNORE_TEST:
      return
    def testSpecies(df, species=None, is_siglvl=False):
      # Ensures that mutations are for the correct species
      if species is not None:
        constraints = [lambda r: r[cn.SPECIES] == species]
      else:
        constraints = []
      self.corr = GenomeCorrelation(constraints=constraints,
          is_test=True, is_siglvl=is_siglvl)
      df = self.corr.makeCorrelationDF()
      if species is None:
        trues = [c[0] in [cn.SPECIES_MIX_DVH, cn.SPECIES_MIX_MMP]
                 for c in df.columns]
      else:
        trues = [c[0] == species for c in df.columns]
      if not all(trues):
        import pdb; pdb.set_trace()
      self.assertTrue(all(trues))
      if is_siglvl:
        df_trues = df.applymap(
            lambda v: ((v <= 1) and (v >= 0) )
                      or np.isclose(v, 1.0))
        trues = [[x for x in c] for c in df_trues.values]
        if not all(trues):
          import pdb; pdb.set_trace()
        self.assertTrue(all(trues))
      # Verify that the correlation  matrix DF is symmetric
      for col1 in df.columns:
        for col2 in df.columns:
          if not np.isclose(df[col1][col2], df[col2][col1]):
            import pdb; pdb.set_trace()
          self.assertTrue(np.isclose(df[col1][col2], df[col2][col1]))

    df = self.corr.makeCorrelationDF()
    self.assertTrue(helpers.isValidDataFrame(df, df.columns))
    #
    for is_siglvl in [False, True]:
      testSpecies(df, is_siglvl=is_siglvl)
      testSpecies(df, species=cn.SPECIES_MIX_DVH, is_siglvl=is_siglvl)
      testSpecies(df, species=cn.SPECIES_MIX_MMP, is_siglvl=is_siglvl)

  def testMakeCorrelationDF2(self):
    if IGNORE_TEST:
      return
    constraints=[
            lambda r: r[cn.SPECIES] == cn.SPECIES_MIX_DVH,
            lambda r: r[cn.LINE] == 'HA2'
    ]
    self.corr = GenomeCorrelation(constraints=constraints,
        is_test=True, is_siglvl=True)
    df = self.corr.makeCorrelationDF()

  def testPlotMutationPairsHist(self):
    if IGNORE_TEST:
      return
    df = self.corr.plotMutationPairsHist(cn.SPECIES_MIX_DVH)
    self.assertTrue(helpers.isValidDataFrame(df, 
        [gc.MUTE1, gc.MUTE2, cn.COUNT]))

  def testPlotGroubyCountHist(self):
    if IGNORE_TEST:
      return
    df = self.corr.plotGroupbyCountHist(cn.KEY_ISOLATE,
        cn.KEY_MUTATION, cn.SPECIES_MIX_DVH)
    self.assertTrue(helpers.isValidDataFrame(df, 
        [cn.KEY_MUTATION, cn.KEY_ISOLATE]))
    df = self.corr.plotGroupbyCountHist(cn.KEY_MUTATION,
        cn.KEY_ISOLATE, cn.SPECIES_MIX_MMP)
    self.assertTrue(helpers.isValidDataFrame(df, 
        [cn.KEY_MUTATION, cn.KEY_ISOLATE]))

  def testMakeCountDF(self):
    if IGNORE_TEST:
      return
    df = self.corr.makeCountDF()
    self.assertTrue(helpers.isValidDataFrame(df,
        [cn.LINE, cn.SPECIES, gc.COUNT_MUTATION, gc.COUNT_ISOLATE]))

  def testPlotHeatmap(self):
    if IGNORE_TEST:
      return
    # Only smoke test
    def test(is_siglvl):
      self.corr = GenomeCorrelation(
          is_test=True, is_siglvl=is_siglvl)
      df_corr = self.corr.makeCorrelationDF()
      self.corr.plotHeatmap(df_corr)
    #
    for siglvl in [True, False]:
      test(siglvl)

  def testMakeBaseDF(self):
    if IGNORE_TEST:
      return
    df = GenomeCorrelation.makeBaseDF()
    self.assertTrue(helpers.isValidDataFrame(
        df, [
        cn.KEY_MUTATION, cn.GENE_ID, cn.GGENE_ID,
        cn.KEY_ISOLATE, cn.COUNT, cn.GENE_POSITION,
        cn.COMMUNITY, cn.SPECIES, cn.LINE],
        nan_columns=[cn.GENE_ID]))


class TestFunctions(unittest.TestCase):

  def testMakeSiglvlDF(self):
    def test(df, is_symmetric=True):
      # Tests if symmetric
      self.assertTrue(df.equals(df.transpose()) == is_symmetric)
    #
    df_15_15 = genome_correlation.makeSiglvlDF(DF15, DF15)
    test(df_15_15)
    df_152_15 = genome_correlation.makeSiglvlDF(DF152, DF15)
    test(df_152_15, is_symmetric=False)
    df_15_152 = genome_correlation.makeSiglvlDF(DF15, DF152)
    test(df_15_152, is_symmetric=False)
    self.assertTrue(df_15_152.equals(df_152_15.transpose()))
    

if __name__ == '__main__':
    unittest.main()
