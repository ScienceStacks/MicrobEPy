from microbepy.common import constants as cn
from microbepy.common import helpers
from microbepy.common.study_context import StudyContext
from microbepy.data.model_data_provider import ModelDataProvider
from microbepy.plot.mutation_plot import MutationLinePlot

import numpy as np
import pandas as pd
import unittest

IGNORE_TEST = False
IS_PLOT = False
PROVIDER = ModelDataProvider(StudyContext(depvar=cn.RATE,
    mutation_column=cn.GGENE_ID))
PROVIDER.do()
SMALL = 1e-8


class TestMutationLinePlot(unittest.TestCase):

  def init(self):
    self.mutation_plot = MutationLinePlot(is_plot=IS_PLOT)

  def setUp(self):
    if IGNORE_TEST:
      return
    self.init()

  def testPlotLine(self):
    # Smoke tests
    if IGNORE_TEST:
      return
    self.mutation_plot.plotLine(cn.TRANSFER_DEFAULT)

  def testPlotTransfers(self):
    if IGNORE_TEST:
      return
    df = self.mutation_plot.plotTransfers()
    self.mutation_plot.plotTransfers(is_cluster_mutations=False)
    self.assertTrue(helpers.isValidDataFrame(df, df.columns))
    self.assertTrue(cn.TRANSFER in df.columns)
    _ = self.mutation_plot.plotTransfers(is_unit_fraction=True)

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

  def testPlotCofraction(self):
    # Smoke test
    if IGNORE_TEST:
      return
    self.mutation_plot.plotCofraction(transfer=cn.TRANSFER_DEFAULT,
        other_transfer=15, is_center_colorbar=False,
        is_differenced=True, is_compress=True)
    self.mutation_plot.plotCofraction(is_center_colorbar=False)
    self.mutation_plot.plotCofraction(transfer=cn.TRANSFER_DEFAULT,
        other_transfer=15, is_center_colorbar=False,
        is_differenced=True)
    self.mutation_plot.plotCofraction(transfer=cn.TRANSFER_DEFAULT,
        other_transfer=15, is_center_colorbar=False,
        is_differenced=False)

  def testPlotCofractions(self):
    # Smoke test
    if IGNORE_TEST:
      return
    self.init()
    self.mutation_plot.plotCofractions(threshold_frac=0.3,
        is_time_lag=True, is_differenced=True)
    self.mutation_plot.plotCofractions(threshold_frac=0.3)



if __name__ == '__main__':
    unittest.main()
