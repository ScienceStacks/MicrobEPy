from microbepy.common import constants as cn
from microbepy.common import helpers
from microbepy.common.study_context import StudyContext
from microbepy.data.model_data_provider import ModelDataProvider
from microbepy.plot.mutation_plot  \
    import MutationIsolatePlot, MutationLinePlot

import numpy as np
import pandas as pd
import unittest

IGNORE_TEST = False
IS_PLOT = False
PROVIDER = ModelDataProvider(StudyContext(depvar=cn.RATE,
    mutation_column=cn.GGENE_ID))
PROVIDER.do()
SMALL = 1e-8


class TestMutationIsolatePlot(unittest.TestCase):

  def setUp(self):
    self.mutation_plot = MutationIsolatePlot(provider=PROVIDER,
        is_plot=IS_PLOT)

  def testMakeMutationIsolateDF(self):
    if IGNORE_TEST:
      return
    df = self.mutation_plot._makeMutationIsolateDF(cn.SPECIES_MIX_DVH)
    self.assertTrue(helpers.isValidDataFrame(df, df.columns))

  def testPlot(self):
    if IGNORE_TEST:
      return
    # Smoke tests
    self.mutation_plot.plot(cn.SPECIES_MIX_DVH)


class TestMutationIsolatePlot(unittest.TestCase):

  def setUp(self):
    self.mutation_plot = MutationLinePlot(is_plot=IS_PLOT)

  def testMakeLineDF(self):
    if IGNORE_TEST:
      return
    df_DVH = self.mutation_plot._makeLineDF(
          species=cn.SPECIES_MIX_DVH, 
          transfer=cn.TRANSFER_DEFAULT)
    df_MMP = self.mutation_plot._makeLineDF(
          species=cn.SPECIES_MIX_MMP, 
          transfer=cn.TRANSFER_DEFAULT)
    df_both = self.mutation_plot._makeLineDF(
          species=None,
          transfer=cn.TRANSFER_DEFAULT)
    self.assertEqual(len(df_DVH) + len(df_MMP), len(df_both))
    for transfer in self.mutation_plot.getTransfers():
      df = self.mutation_plot._makeLineDF(
          species=cn.SPECIES_MIX_DVH, 
          transfer=transfer)
      self.assertTrue(helpers.isValidDataFrame(df, df.columns))

  def testPlotLine(self):
    # Smoke tests
    if IGNORE_TEST:
      return
    self.mutation_plot.plotLine(cn.SPECIES_MIX_DVH, 
        cn.TRANSFER_DEFAULT)

  def testGetFrequentMutations(self):
    if IGNORE_TEST:
      return
    mutations_2 = self.mutation_plot._getFrequentMutations(min_lines=2)
    mutations_3 = self.mutation_plot._getFrequentMutations(min_lines=3)
    self.assertGreater(len(mutations_2), len(mutations_3))
    #
    mutations_dvh = self.mutation_plot._getFrequentMutations(
        species=cn.SPECIES_MIX_DVH)
    mutations_mmp = self.mutation_plot._getFrequentMutations(
        species=cn.SPECIES_MIX_MMP)
    self.assertGreater(len(mutations_dvh), len(mutations_mmp))

  def testGetLines(self):
    if IGNORE_TEST:
      return
    lines = self.mutation_plot.getLines()
    self.assertGreater(len(lines), 0)
    #
    lines_mmp = self.mutation_plot.getLines(cn.SPECIES_MIX_MMP)
    lines_dvh = self.mutation_plot.getLines(cn.SPECIES_MIX_DVH)
    self.assertEqual(set(lines), set(lines_mmp).union(lines_dvh))

  def testPlotTransfers(self):
    if IGNORE_TEST:
      return
    df = self.mutation_plot.plotTransfers()
    self.mutation_plot.plotTransfers(is_cluster_mutations=False)
    self.assertTrue(helpers.isValidDataFrame(df, df.columns))
    self.assertTrue(cn.TRANSFER in df.columns)
    _ = self.mutation_plot.plotTransfers(is_unit_fraction=True)

  def testOrderMutations(self):
    if IGNORE_TEST:
      return
    mutations_all = self.mutation_plot._orderMutations()
    mutations_DVH = self.mutation_plot._orderMutations(
        species=cn.SPECIES_MIX_DVH)
    mutations_MMP = self.mutation_plot._orderMutations(
        species=cn.SPECIES_MIX_MMP)
    mutations_union = set(mutations_DVH).union(mutations_MMP)
    empty = mutations_union.symmetric_difference(mutations_all)
    self.assertEqual(len(empty), 0)

  def testMakeMutationSiglvlMatrix(self):
    if IGNORE_TEST:
      return
    df_matrix = self.mutation_plot._makeMutationSiglvlMatrix()
    self.assertEqual(set(df_matrix.columns), set(df_matrix.index))
    values = []
    for col in df_matrix:
      self.assertTrue(all([(v <= 1 + SMALL) and (v >= -SMALL) 
          for v in df_matrix[col]]))

  def testPlotSiglvlDF(self):
    if IGNORE_TEST:
      return
    df = self.mutation_plot._plotSiglvlDF(max_siglvl=0.5)
    self.assertTrue(helpers.isValidDataFrame(df, df.columns,
        nan_columns=df.columns))
    df = self.mutation_plot._plotSiglvlDF(transfer=15,
        other_transfer=45,
        max_siglvl=0.5)
    self.assertTrue(helpers.isValidDataFrame(df, df.columns,
        nan_columns=df.columns))

  def testPlotSiglvl(self):
    if IGNORE_TEST:
      return
    # Smoke test
    self.mutation_plot.plotSiglvl(transfer=15,
        other_transfer=45, max_siglvl=0.05)
    self.mutation_plot.plotSiglvl(transfer=45,
        other_transfer=76, max_siglvl=0.05, is_center_colorbar=False)
    self.mutation_plot.plotSiglvl(transfer=76,
        other_transfer=118, max_siglvl=0.05)
    self.mutation_plot.plotSiglvl(transfer=118,
        other_transfer=152, max_siglvl=0.05)
    self.mutation_plot.plotSiglvl(max_siglvl=0.05)

  def testPlotSiglvl(self):
    if IGNORE_TEST:
      return
    # Smoke test
    self.mutation_plot.plotSiglvls(max_siglvl=0.05)
    self.mutation_plot.plotSiglvls(max_siglvl=0.05, is_time_lag=True)


if __name__ == '__main__':
    unittest.main()
