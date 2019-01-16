from microbepy.common import helpers
from microbepy.common import util
from microbepy.common import constants as cn
from microbepy.plot import regression_plot as rp
from microbepy.plot.util_plot import PlotParms
from microbepy.model import isolate_regression as ir

import numpy as np
import pandas as pd
import unittest

IGNORE_TEST = False
IS_TEST = True
DEPVAR = cn.YIELD


################ CLASS UNDER TEST ############
class TestFunctions(unittest.TestCase):

  def setUp(self):
    self.df = ir.NonParametricIsolateRegression.makeEstimateDFS(
        depvar=DEPVAR)[cn.RESIDUAL]

  def testPlotFilteredEstimatedObservedLeaveOutIsolates(self):
    if IGNORE_TEST:
      return
    self.rplot = rp.RegressionPlot(self.df, is_test=IS_TEST)
    depvar = cn.RATE
    line = 'HA2'
    isolate_pair = ("HA2.152.08.03.D.CI", "HA2.152.08.03.M.CI")
    self.rplot.plotFilteredEstimatedObserved(
        ir.AncestralPairingIsolateRegression, depvar, max_std=3, 
            line=line, leave_out=cn.KEY_CULTURE,
            leave_out_isolates=[isolate_pair])

  def testPlotFilteredEstimatedObservedCulture(self):
    if IGNORE_TEST:
      return
    self.rplot = rp.RegressionPlot(self.df, is_test=IS_TEST)
    depvar = cn.RATE
    line = 'HA2'
    self.rplot.plotFilteredEstimatedObserved(
        ir.AncestralPairingIsolateRegression, depvar, max_std=3, 
            line=line, leave_out=cn.KEY_CULTURE)

  def testPlotFilteredEstimatedObservedIsolate(self):
    if IGNORE_TEST:
      return
    self.rplot = rp.RegressionPlot(self.df, is_test=IS_TEST)
    depvar = cn.RATE
    line = 'HR2'
    self.rplot.plotFilteredEstimatedObserved(
        ir.AncestralPairingIsolateRegression, depvar, max_std=3,
            line=line, leave_out=cn.KEY_ISOLATE)

  def testPlotFilteredEstimatedObserved2(self):
    if IGNORE_TEST:
      return
    self.rplot = rp.RegressionPlot(self.df, is_test=IS_TEST)
    depvar = cn.RATE
    self.rplot.plotFilteredEstimatedObserved(
        ir.NonParametricIsolateRegression, depvar, max_std=3)
    

if __name__ == '__main__':
    unittest.main()
