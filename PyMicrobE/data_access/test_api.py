
import __init__
import util
from api import Api
import constants as cn
import api as aa
import constants as cn
import helpers

import numpy as np
import os
import pandas as pd
import unittest
import random, string

IGNORE_TEST = False


########################################
class TestAPI(unittest.TestCase):


  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.api = Api()
    trues = [isinstance(df, pd.DataFrame) 
             for df in list(self.api.dfs.values())]
    self.assertTrue(all(trues))

  def testReadSQL(self):
    if IGNORE_TEST:
      return
    self.api = Api()
    for table in [cn.TABLE_MUTATION, cn.TABLE_CULTURE, cn.TABLE_ISOLATE]:
      schema = cn.TABLE_SCHEMAS.getSchema(table)
      sql_cmd = "SELECT * from %s;" % table
      df = self.api.readSQL(sql_cmd)
      difference = set(schema.columns).symmetric_difference(
          df.columns)
      self.assertEqual(len(difference), 0)
      self.assertEqual(len(df.index),
          len(self.api.dfs[table].index))

  def testMakeDF(self):
    if IGNORE_TEST:
      return
    self.api = Api()
    def testSingleTable(table_name):
      schema = self.api.schemas[table_name]
      if schema.name == cn.TABLE_CULTURE_ISOLATE_MUTATION:
        no_null_columns = []
      else:
        no_null_columns = schema.key
      df_raw = util.cleanDF(self.api.dfs[table_name])
      df = self.api.makeDF(schema.key, 
          no_null_columns=no_null_columns)
      self.assertEqual(len(df.index), 
          len(self.api.dfs[table_name].index))
    #
    def testColumnCombinations(columns, expected_length):
      df = self.api.makeDF(columns)
      self.assertEqual(len(df.index), expected_length)
    #
    df_mute_gene = self.api.dfs[cn.TABLE_MUTATION][
        [cn.GENE_ID, cn.CHANGED_NT]]
    sel = pd.Series([not util.isNull(x) 
        for x in df_mute_gene[cn.GENE_ID]])
    df_mute_gene = df_mute_gene.loc[sel]
    df_mute_gene.drop_duplicates(inplace=True)
    testColumnCombinations([cn.GENE_ID, cn.CHANGED_NT], 
        len(df_mute_gene.index))
    # Verify that individual tables are handled correctly
    test_tables  = [cn.TABLE_MUTATION, cn.TABLE_ISOLATE,
        cn.TABLE_CULTURE, cn.TABLE_ISOLATE_MUTATION_LINK,
        cn.TABLE_CULTURE_ISOLATE_LINK]
    for table_name in test_tables:
      testSingleTable(table_name)
    #
    testColumnCombinations([cn.KEY_MUTATION, cn.KEY_ISOLATE],
        len(self.api.dfs[cn.TABLE_ISOLATE_MUTATION_LINK].index))
    testColumnCombinations([cn.KEY_CULTURE, cn.KEY_ISOLATE],
        len(self.api.dfs[cn.TABLE_CULTURE_ISOLATE_LINK].index))
    testColumnCombinations([cn.GENE_ID, cn.KEY_MUTATION],
        len(self.api.dfs[cn.TABLE_MUTATION].index))

  def testMakeDFOptions(self):
    if IGNORE_TEST:
      return
    self.api = Api()
    length = len(self.api.dfs[cn.TABLE_MUTATION].index)
    # No options
    columns = [cn.KEY_MUTATION, cn.MUTATION]
    df = self.api.makeDF(columns)
    self.assertEqual(len(df.index), length)
    self.assertEqual(set(df.columns), set(columns))
    # Species mix and effects
    columns = [cn.KEY_MUTATION, cn.MUTATION]
    df = self.api.makeDF(columns,
        effects=cn.EFFECTS, species_mixes=cn.SPECIES_MIXES)
    self.assertGreater(len(df.index), 0)
    self.assertEqual(set(df.columns), set(columns))
    # Constraint without aux_columns
    sel_key = self.api.dfs[
       cn.TABLE_MUTATION][cn.KEY_MUTATION].tolist()[0]
    constraint = lambda r: r[cn.KEY_MUTATION] == sel_key
    df = self.api.makeDF(columns, constraints=[constraint])
    self.assertEqual(len(df.index), 1)
    # Constraint with an aux_column
    constraint = lambda r: r[cn.EFFECT] == "DUMMY"
    aux_columns = [cn.EFFECT]
    df = self.api.makeDF(columns, constraints=[constraint],
        aux_columns=aux_columns)
    self.assertEqual(len(df.index), 0)

  def testConsistency1(self):
    """
    Verify that correctly get cultures that should have species
    pairs.
    """
    if IGNORE_TEST:
      return
    self.api = Api()
    columns = [cn.KEY_CULTURE, cn.KEY_ISOLATE]
    df_both = self.api.makeDF(columns=columns,
        aux_columns = [cn.SPECIES_MIX],
        species_mixes=cn.SPECIES_MIX_BOTH)
    df_dvh = self.api.makeDF(columns=columns, 
        species_mixes=cn.SPECIES_MIX_BOTH,
        aux_columns = [cn.SPECIES, cn.SPECIES_MIX],
        constraints=[lambda r: r[cn.SPECIES]==cn.SPECIES_MIX_DVH])
    df_mmp = self.api.makeDF(
        columns=columns,
        aux_columns=[cn.SPECIES, cn.SPECIES_MIX],
        species_mixes=cn.SPECIES_MIX_BOTH,
        constraints=[lambda r: r[cn.SPECIES]==cn.SPECIES_MIX_MMP])
    self.assertEqual(2*len(df_dvh.index), len(df_both.index))
    self.assertEqual(2*len(df_mmp.index), len(df_both.index))

  def testConsistency2(self):
    """
    Verify that correctly get cultures that should have species
    pairs.
    """
    KEY_ISOLATE_DVH = "key_isolate_dvh"
    KEY_ISOLATE_MMP = "key_isolate_mmp"
    self.api = Api()
    columns = [cn.KEY_CULTURE, cn.RATE, cn.YIELD, cn.KEY_ISOLATE, 
        cn.LINE, cn.SPECIES]
    df_dvh = self.api.makeDF(columns=columns, 
        effects=cn.EFFECTS,
        species_mixes=cn.SPECIES_MIX_BOTH,
        constraints=[lambda r: r[cn.SPECIES]==cn.SPECIES_MIX_DVH])
    df_dvh.rename(columns={
        cn.KEY_ISOLATE: KEY_ISOLATE_DVH}, inplace=True)
    df_mmp = self.api.makeDF(
        columns=[cn.KEY_CULTURE, cn.KEY_ISOLATE], 
        aux_columns=[cn.SPECIES],
        effects=cn.EFFECTS,
        species_mixes=cn.SPECIES_MIX_BOTH,
        constraints=[lambda r: r[cn.SPECIES]==cn.SPECIES_MIX_MMP])
    df_mmp.rename(columns={
      cn.KEY_ISOLATE: KEY_ISOLATE_MMP}, inplace=True)
    df = df_dvh.merge(df_mmp, on=cn.KEY_CULTURE, how='inner')
    self.assertEqual(len(df_dvh.index), len(df.index))
    self.assertEqual(len(df_mmp.index), len(df.index))
    
    
if __name__ == '__main__':
    unittest.main()
