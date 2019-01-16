from microbepy.common import constants as cn
from microbepy.common.isolate import Isolate
from microbepy.common import isolate
from microbepy.common import util

import numpy as np
import os
import pandas as pd
import unittest
import random, string

IGNORE_TEST = False
LINE = 'HR'
LINE_REPLICA_NUMBER = "2"
LINE_REPLICA = LINE + LINE_REPLICA_NUMBER
EPD = LINE_REPLICA + ".152.02"
ISOLATE_NUMBER = "01"
DVH_ISOLATE = EPD + "." + cn.SPECIES_MIX_DVH + ISOLATE_NUMBER
MMP_ISOLATE = EPD + "." + cn.SPECIES_MIX_MMP + ISOLATE_NUMBER
WT_DVH_ISOLATE = "WT.D01"
WT_MMP_ISOLATE = "WT.M01"


########################################
class TestIsolate(unittest.TestCase):

  def setUp(self):
    if IGNORE_TEST:
      return
    pass

  def testGetSpecies(self):
    if IGNORE_TEST:
      return
    def test(isolate_stg, species):
      isolate = Isolate.create(isolate_stg)
      self.assertEqual(isolate.species, species)

    test(DVH_ISOLATE, cn.SPECIES_MIX_DVH)
    test(MMP_ISOLATE, cn.SPECIES_MIX_MMP)
    test(WT_DVH_ISOLATE, cn.SPECIES_MIX_DVH)
    test(WT_MMP_ISOLATE, cn.SPECIES_MIX_MMP)

  def testIsWTorAN(self):
    if IGNORE_TEST:
      return
    self.assertFalse(Isolate.isAN(DVH_ISOLATE))
    with self.assertRaises(ValueError):
      self.assertFalse(Isolate.isAN("garbage"))
    self.assertTrue(Isolate.isAN("AN.D01"))

  def testValidIsolate(self):
    # Tests that particular isolate constructions are valid.
    # These are all smoke tests. There should be no exceptions.
    Isolate(line=cn.LINE_ANC)
    for line in [cn.LINE_AN1, cn.LINE_AN2]:
      for species in [cn.SPECIES_MIX_DVH, cn.SPECIES_MIX_MMP]:
        Isolate(line=line, species=species)

  def testGetLine(self):
    if IGNORE_TEST:
      return
    self.assertEqual(Isolate.create(DVH_ISOLATE).line[0:2], LINE)

  def testIsSpecies(self):
    if IGNORE_TEST:
      return
    dvh_isolate = 'UA2.152.01.D01'
    mmp_isolate = 'UA2.152.01.M01'
    non_isolate = 'UA2.152.01.X01'
    self.assertTrue(Isolate.isSpecies(dvh_isolate, 
        cn.SPECIES_MIX_DVH))
    self.assertTrue(Isolate.isSpecies(mmp_isolate, 
        cn.SPECIES_MIX_MMP))
    with self.assertRaises(ValueError):
      Isolate.isSpecies(non_isolate, cn.SPECIES_MIX_MMP)

  def testIsEPD(self):
    if IGNORE_TEST:
      return
    self.assertTrue(Isolate.isEPD('UA2.152.01.*'))

  def testIsLine(self):
    if IGNORE_TEST:
      return
    isolate_stg = 'UA2.152.01.01.M.SC'
    self.assertTrue(Isolate.isLine(isolate_stg, 'UA'))
    self.assertTrue(Isolate.isLine(isolate_stg, 'UA2'))
    isolate_stg = 'AN.*.*.01.M.*'
    self.assertTrue(Isolate.isLine(isolate_stg, 'AN'))

  def testGetEPDCommunity(self):
    if IGNORE_TEST:
      return
    def test(isolate_stg, expected):
      isolate = Isolate.create(isolate_stg)
      self.assertEqual(isolate.getEPDCommunity(),
        expected)
    #
    test('UA2.152.01.*.*.*', 'UA2.152.01')
    test('UA2.152.01.03.D.CI', 'UA2.152.01')
    test('UA2.152.11.03.D.CI', 'UA2.152.11')

  def testGetClonePairingID(self):
    if IGNORE_TEST:
      return
    def test(isolate_stg, expected):
      isolate = Isolate.create(isolate_stg)
      self.assertEqual(isolate.getClonePairingID(),
        expected)
    #
    test('UA2.152.01.01.*.*', 'UA2.152.01.01')
    test('UA2.152.01.03.D.CI', 'UA2.152.01.03')
    test('UA2.152.11.03.D.CI', 'UA2.152.11.03')

  def testGetCommunity(self):
    if IGNORE_TEST:
      return
    def test(isolate_stg, expected):
      isolate = Isolate.create(isolate_stg)
      self.assertEqual(isolate.getCommunity(),
        expected)
    #
    test('UA2.152.01.*.*.*', cn.ISOLATE_EPD)
    test('UA2.152.01.01.M.CI', cn.ISOLATE_CLONE)
    test('UA2.152.01.01.M.SC', cn.ISOLATE_CLONE)
    test('UA2.152.*.*.*.*', cn.ISOLATE_LINE)
    test('AN.*.*.01.D.*', cn.ISOLATE_LINE)
    test('WT.*.*.01.D.*', cn.ISOLATE_LINE)
    test('*.*.*.*.*.*', cn.ISOLATE_UNKNOWN)
    test(np.nan, cn.ISOLATE_UNKNOWN)

  def testCreate(self):
    if IGNORE_TEST:
      return
    good_isolate_stgs = [
        'UA2.152.01.*.*.*',
        'UA2.152.01.01.M.CI',
        'UA2.152.*.*.*.*',
        'AN.*.*.01.D.*',
        '*.*.*.*.*.*',
        'WT.*.*.01.M.*']
    for stg in good_isolate_stgs:
      isolate = Isolate.create(stg)
      self.assertEqual(str(isolate), stg)
    bad_isolate_stgs = [
        'AN.*.*.*.*.*',
        'WT.*.*.*.*.*',
        'UA2.*.01.01.D.*',
        ]
    for stg in bad_isolate_stgs:
      try:
        _ = Isolate.create(stg)
        import pdb; pdb.set_trace()
      except:
        pass
      with self.assertRaises(ValueError):
        _ = Isolate.create(stg)

  def testCheckDefaults(self):
    if IGNORE_TEST:
      return
    def setValue(v):
      if v is None:
        return []
      else: 
        return [v]
    #
    def test(default_values, non_default_values, expected):
      """
      :param list default_values:
      :param list non_default_values:
      :param bool expected:
      """
      dv_list = []
      for dv in default_values:
        dv_list.extend(setValue(dv))
        ndv_list = []
        for ndv in non_default_values:
          ndv_list.extend(setValue(ndv))
          self.assertEqual(isolate.checkDefaults(dv_list, ndv_list),
              expected)
    #
    default_values = [None, np.nan, cn.ISOLATE_DEFAULT]
    non_default_values = [None, cn.LINE, cn.EPD_ID]
    test(default_values, non_default_values, True)
    default_values = [np.nan, cn.ISOLATE_DEFAULT]
    non_default_values = [cn.LINE, cn.EPD_ID]
    test(non_default_values, default_values, False)
    #
    self.assertFalse(isolate.checkDefaults(
        ['01', '*', '*', '*'], ['UA2', '152']))

  def testLine(self):
    if IGNORE_TEST:
      return
    isolate = "UE3.152.09.*.*.*"
    line = Isolate.create(isolate).line
    self.assertEqual(line, cn.LINE_UE3)
      
    

if __name__ == '__main__':
    unittest.main()
