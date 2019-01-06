
import __init__
import constants as cn
import regression_plot as rp
from util_plot import PlotParms

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import unittest

IGNORE_TEST = False
IS_PLOT = False
ATTR = "dummy"


############ TEST CLASS ###############
class PlotTester(object):

  def __init__(self):
    self.xvals = range(10)
    self.yvals = range(10)

  def do(self, parms):
    plt.plot(self.xvals, self.yvals)
    parms.do(is_plot=IS_PLOT)


################ CLASS UNDER TEST ############
class TestPlotParms(unittest.TestCase):

  def setUp(self):
    self.plot = PlotTester()

  def testPlotParms(self):
    if IGNORE_TEST:
      return
    parms = PlotParms()
    parms[cn.PLT_TITLE] = "title"
    parms[cn.PLT_XLIM] = [0, 5]
    parms[cn.PLT_YLIM] = [0, 5]
    self.plot.do(parms)

  def testIsTrue(self):
    if IGNORE_TEST:
      return
    parms = PlotParms()
    self.assertFalse(parms.isTrue('a'))
    parms['a'] = False
    self.assertFalse(parms.isTrue('a'))
    parms['a'] = True
    self.assertTrue(parms.isTrue('a'))

  def testSetTrueIfAbsent(self):
    if IGNORE_TEST:
      return
    parms = PlotParms()
    self.assertFalse(parms.isTrue(ATTR))
    parms.setTrueIfAbsent(ATTR)
    self.assertTrue(parms.isTrue(ATTR))

if __name__ == '__main__':
    unittest.main()
