import __init__
import constants as cn
from data_model import DataModel
import data_model as dm
import helpers
from isolate import Isolate
import util
import util_data_access

import numpy as np
import os
import pandas as pd
import unittest

IGNORE_TEST = False

class TestDataModel(unittest.TestCase):

  def setUp(self):
    self.data_model = DataModel(output_directory=cn.TEST_PATH)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(
        isinstance(self.data_model._dvh_mutation_file, str))
    self.assertIsNotNone(self.data_model._dvh_mutation_file)

  def testMakePairingDF(self):
    if IGNORE_TEST:
      return
    df = self.data_model._makePairingDF()
    self.assertTrue(helpers.isValidDataFrame(df,
        expected_columns=[cn.SPECIES, cn.KEY_ISOLATE, cn.SAMPLE],
        key=cn.KEY_ISOLATE, min_rows=70))

  def testMakeMutationDF(self):
    if IGNORE_TEST:
      return
    self.data_model._makeCombinationDFs()
    df = self.data_model._makeMutationDF()
    expecteds = cn.TABLE_SCHEMAS.getSchema(cn.TABLE_MUTATION).columns
    key = cn.TABLE_SCHEMAS.getSchema(cn.TABLE_MUTATION).key
    self.assertTrue(helpers.isValidDataFrame(df,
        expected_columns=expecteds,
        key=key,
        nan_columns=expecteds,
        ))
    self.assertGreater(
        sum([1 if x else 0 for x in df[cn.IS_AN_MUTATION]]), 0)

  def testMakeIsolatesForMutations(self):
    if IGNORE_TEST:
      return
    def test(mutation_file, species):
      df_mutation = util.makeDataframeFromXlsx(
          util.getSequenceDataPath(mutation_file))
      df_mutation[cn.SPECIES] = species
      isolate_stgs = self.data_model._makeIsolatesForMutations(
          df_mutation)
      for stg in isolate_stgs:
        try:
          _ = Isolate.create(stg)
        except:
          # Bad isolate string
          self.assertTrue(False)
      self.assertEqual(len(isolate_stgs), len(df_mutation.index))

    test(dm.DVH_MUTATION_FILE, cn.SPECIES_MIX_DVH)
    test(dm.MMP_MUTATION_FILE, cn.SPECIES_MIX_MMP)

  def testMakeCultureDF(self):
    if IGNORE_TEST:
      return
    self.data_model._makeCombinationDFs()
    df_culture_isolate = self.data_model._makeCultureIsolateDF()
    df = self.data_model._makeCultureDF()
    expecteds = cn.TABLE_SCHEMAS.getSchema(cn.TABLE_CULTURE).columns
    key = cn.TABLE_SCHEMAS.getSchema(cn.TABLE_CULTURE).key
    self.assertTrue(helpers.isValidDataFrame(df,
        expected_columns=expecteds,
        nan_columns = expecteds,
        key=key))

  def testMakeSpeciesMixInCulture(self):
    if IGNORE_TEST:
      return
    df_culture_isolate = self.data_model._makeCultureIsolateDF()
    df_mix = dm.makeSpeciesMixInCulture(df_culture_isolate)
    self.assertTrue(helpers.isValidDataFrame(df_mix,
        expected_columns=[cn.KEY_CULTURE, cn.SPECIES_MIX],
        key=cn.KEY_CULTURE))
    self.assertTrue(
        set(df_mix[cn.SPECIES_MIX].unique()).issubset(cn.SPECIES_MIXES))

  def testMakeIsolateMutationDF(self):
    if IGNORE_TEST:
      return
    df = self.data_model._makeIsolateMutationDF()
    expecteds = list(cn.TABLE_SCHEMAS.getSchema(
        cn.TABLE_MUTATION).columns)
    expecteds.remove(cn.IS_AN_MUTATION)
    expecteds.remove(cn.GENE_POSITION)
    expecteds.extend([cn.KEY_ISOLATE])
    nan_columns=[cn.INITIAL_NT, cn.CHANGED_NT, cn.AA_CHANGE,
        cn.GENE_ID, cn.CODON_CHANGE, 'region']
    self.assertTrue(
        helpers.isValidDataFrame(df, expected_columns=expecteds,
        nan_columns=nan_columns))

  def testMakeIsolateDF(self):
    if IGNORE_TEST:
      return
    self.data_model._makeCombinationDFs()
    df_isolate_mutation = self.data_model._makeIsolateMutationDF()
    df_culture_isolate = self.data_model._makeCultureIsolateDF()
    df = self.data_model._makeIsolateDF()
    expecteds = cn.TABLE_SCHEMAS.getSchema(cn.TABLE_ISOLATE).columns
    key = cn.TABLE_SCHEMAS.getSchema(cn.TABLE_ISOLATE).key
    self.assertTrue(helpers.isValidDataFrame(
        df, 
        expected_columns=expecteds,
        nan_columns=expecteds,
        key=key))
    count = sum([1 if (x and not util.isNull(x)) else 0 
                 for x in df[cn.IS_LOW_COVERAGE_ISOLATE]])
    self.assertEqual(count, len(cn.LOW_COVERAGE_ISOLATES))

  def testMakeIsolateMutationLink(self):
    if IGNORE_TEST:
      return
    self.data_model._makeCombinationDFs()
    df_culture_isolate = self.data_model._makeCultureIsolateDF()
    df_isolate_mutation = self.data_model._makeIsolateMutationDF()
    df = self.data_model._makeIsolateMutationLinkDF()
    valid_dict = {
        cn.KEY_ISOLATE: lambda x: util.isStr(x),
        cn.KEY_MUTATION: lambda x: util.isNull(x) or util.isStr(x),
        cn.FREQ: lambda x: isinstance(x, float),
        cn.GATK: lambda x: isinstance(x, float),
        cn.SAMTOOLS: lambda x: isinstance(x, float),
        cn.VARSCAN: lambda x: isinstance(x, float),
        }
    expecteds = cn.TABLE_SCHEMAS.getSchema(
        cn.TABLE_ISOLATE_MUTATION_LINK).columns
    key = cn.TABLE_SCHEMAS.getSchema(
        cn.TABLE_ISOLATE_MUTATION_LINK).key
    nan_columns = list(cn.TABLE_SCHEMAS.getSchema(
        cn.TABLE_ISOLATE_MUTATION_LINK).columns)
    nan_columns.remove(cn.KEY_ISOLATE)
    self.assertTrue(helpers.isValidDataFrame(
        df, 
        nan_columns=nan_columns,
        expected_columns=expecteds,
        key=key,
        valid_dict=valid_dict))

  def testMakeCultureIsolateLinkDF(self):
    if IGNORE_TEST:
      return
    self.data_model._makeCombinationDFs()
    df_culture_isolate = self.data_model._makeCultureIsolateDF()
    df_isolate_mutation = self.data_model._makeIsolateMutationDF()
    df = self.data_model._makeCultureIsolateLinkDF()
    expecteds = cn.TABLE_SCHEMAS.getSchema(
        cn.TABLE_CULTURE_ISOLATE_LINK).columns
    key = cn.TABLE_SCHEMAS.getSchema(
         cn.TABLE_CULTURE_ISOLATE_LINK).key
    self.assertTrue(helpers.isValidDataFrame(
        df, 
        nan_columns = [cn.KEY_CULTURE],
        expected_columns=expecteds,
        key=key))

  def testMakeGeneDescriptionDF(self):
    if IGNORE_TEST:
      return
    self.data_model._makeCombinationDFs()
    df = self.data_model._makeGeneDescriptionDF()
    expecteds = cn.TABLE_SCHEMAS.getSchema(
        cn.TABLE_GENE_DESCRIPTION).columns
    nan_columns = df.columns.tolist()
    self.assertTrue(helpers.isValidDataFrame(
        df,
        nan_columns = nan_columns, 
        expected_columns=expecteds,
        ))

  def testDo(self):
    if IGNORE_TEST:
      return
    # Eliminate any existing output
    for schema in cn.TABLE_SCHEMAS.getSchemas():
      csv_file = "%s.csv" % schema.name
      path = os.path.join(cn.TEST_PATH, csv_file)
      if os.path.isfile(path):
        os.remove(path)
    # Create new output
    self.data_model.do()
    # Validate the new output
    for schema in cn.TABLE_SCHEMAS.getSchemas():
      csv_file = "%s.csv" % schema.name
      path = os.path.join(cn.TEST_PATH, csv_file)
      df = pd.read_csv(path)
      # Check which is the needed key
      columns = df.columns.tolist()
      self.assertTrue(helpers.isValidDataFrame(df,
          expected_columns=schema.columns,
          key=schema.key,
          nan_columns=schema.columns))

  def testMakeCultureIsolateMutationDF(self):
    if IGNORE_TEST:
      return
    def test(df1, df2, schema):
      """
      Checks the equivalence of two dataframes.
      :param pd.DataFrame df1:
      :param pd.DataFrame df2:
      """
      # Check the columns
      columns = set(df1.columns)
      differences = columns.symmetric_difference(df2.columns)
      if not len(differences) == 0:
          import pdb; pdb.set_trace()
      self.assertEqual(len(differences), 0)    
      # Find if the columns differ
      df1 = df1.sort_values(schema.key)
      df2 = df2.sort_values(schema.key)
      rows1 = [r for _,r in df1.iterrows()]
      rows2 = [r for _,r in df2.iterrows()]
      length = min(len(rows1), len(rows2))
      for idx in range(length):
        row1 = rows1[idx]
        row2 = rows2[idx]
        for column in columns:
          val1 = row1[column]
          val2 = row2[column]
          try:
            if util.isStr(val1) or util.isStr(val2):
              val1 = str(val1)
              val2 = str(val2)
          except:
            import pdb; pdb.set_trace()
          if util.isNull(val1) or util.isNull(val2):
            b = util.isNull(val1) and util.isNull(val2)
          elif util.isNumber(val1) or util.isNumber(val2):
            b = np.isclose(float(val1), float(val2))
          else:
            b = val1 == val2
          if not b:
            # Try running data_model.py to generate the csv files
            # and then rerun the tests.
            import pdb; pdb.set_trace()
          self.assertTrue(b)
      # Check rows
      if len(df1.index) != len(df2.index):
        import pdb; pdb.set_trace()
      self.assertEqual(len(df1.index), len(df2.index))

    # Construct the inputs for the universal DF
    self.data_model._makeCombinationDFs()
    self.data_model._makeBaseDataModelDFs()
    # Test DataFrame basics
    df_cim = self.data_model._makeCultureIsolateMutationDF()
    schema = cn.TABLE_SCHEMAS.getSchema(
        cn.TABLE_CULTURE_ISOLATE_MUTATION)
    df_key = df_cim[[cn.KEY_CULTURE, cn.KEY_ISOLATE,
        cn.KEY_MUTATION, cn.INITIAL_NT]]
    dfg = df_key.groupby([cn.KEY_CULTURE, cn.KEY_ISOLATE,
        cn.KEY_MUTATION]).count()
    self.assertTrue(helpers.isValidDataFrame(df_cim,
        expected_columns=schema.columns,
        key = schema.key,
        nan_columns = schema.columns))
    # Verify each table
    tables = [cn.TABLE_MUTATION, cn.TABLE_ISOLATE, cn.TABLE_CULTURE,
        cn.TABLE_CULTURE_ISOLATE_LINK, cn.TABLE_ISOLATE_MUTATION_LINK]
    for table in tables:
      schema = cn.TABLE_SCHEMAS.getSchema(table)
      df = df_cim.copy(deep=True)
      df = df[schema.columns].copy(deep=True)
      df = util.pruneNullRows(df)
      # Compound keys cannot have a nan value
      for key in schema.key:
        sel = [not util.isNull(x) for x in df[key]]
        df = df.loc[sel,:]
      df.drop_duplicates(inplace=True)
      try:
        df1 = util_data_access.readDataModelCSV(schema)
        test(df1, df, schema)
      except FileNotFoundError:
        # Ignore this warning if have *.db file
        print ("**Warning: No data model CSV for %s" % schema.name)
        

  def _setupDF(self):
    # Sets up internal state for data model combination DataFrames
    self.data_model._makeCombinationDFs()
    self.data_model._makeBaseDataModelDFs()
    self.data_model._dfs[cn.TABLE_CULTURE_ISOLATE_MUTATION] =  \
        self.data_model._makeCultureIsolateMutationDF()

  def testMakeGenotypePhenotypeDF(self):
    if IGNORE_TEST:
      return
    self._setupDF()
    df =  self.data_model._makeGenotypePhenotypeDF()
    #
    def test(table, key):
      # Ensure that the length of the constructed
      # dataframe is what's expected compared with the composing
      # tables.
      difference = set(df[key]).difference(
          self.data_model._dfs[table][key])
      self.assertEqual(len(difference), 0)
    #
    columns = cn.TABLE_SCHEMAS.getColumns([cn.TABLE_GENOTYPE_PHENOTYPE])
    self.assertTrue(helpers.isValidDataFrame(df,
        columns, nan_columns=columns,
        key=[cn.KEY_MUTATION, cn.KEY_ISOLATE, cn.KEY_CULTURE]))
    test(cn.TABLE_MUTATION, cn.KEY_MUTATION)
    test(cn.TABLE_CULTURE, cn.KEY_CULTURE)

  def testMakeGenotypeDF(self):
    if IGNORE_TEST:
      return
    self._setupDF()
    df =  self.data_model._makeGenotypeDF()
    #
    def test(table, key):
      # Ensure that the length of the constructed
      # dataframe is what's expected compared with the composing
      # tables.
      difference = set(df[key]).difference(
          self.data_model._dfs[table][key])
      if len(difference) !=  0:
        import pdb; pdb.set_trace()
        pass
    #
    columns = cn.TABLE_SCHEMAS.getColumns([cn.TABLE_GENOTYPE])
    self.assertTrue(helpers.isValidDataFrame(df,
        columns, nan_columns=columns,
        key=[cn.KEY_MUTATION, cn.KEY_ISOLATE]))
    test(cn.TABLE_MUTATION, cn.KEY_MUTATION)
    test(cn.TABLE_ISOLATE_MUTATION_LINK, cn.KEY_ISOLATE)

class TestFunctions(unittest.TestCase):

  def testMakeRawParingDF(self):
    if IGNORE_TEST:
      return
    df = dm.makeRawPairingDF(dm.DVH_PAIRING_FILE)
    self.assertTrue(helpers.isValidDataFrame(df,
        expected_columns=[cn.SAMPLE, dm.KEY_ISOLATE_FRONT, 
        dm.KEY_ISOLATE_BACK]))

  def testMakeKeyMutation(self):
    if IGNORE_TEST:
      return
    def test(mutation_file, species):
      df_mutation = util.makeDataframeFromXlsx(
          util.getSequenceDataPath(mutation_file))
      dm.deleteUnneededMutationColumns(df_mutation)
      df_mutation[cn.SPECIES] = species
      key_mutations = dm.makeKeyMutation(df_mutation)
      self.assertEqual(len(key_mutations), len(df_mutation.index))
      # Ensure that this is a key
      util.deleteColumns(df_mutation, 
          [cn.PREDICTOR, cn.FREQ_PREDICTOR, cn.SAMPLE, 
          cn.VARIANT_ID, cn.EXPERIMENT, cn.FREQ, cn.READ_NUMBER],
          is_drop_duplicates=True)
      self.assertEqual(len(set(key_mutations)), 
          len(df_mutation))
      # Check the elements of the keys
      for stg in key_mutations:
        split = stg.split(cn.MUTATION_SEPARATOR)
        self.assertEqual(len(split), 4)
        self.assertTrue(util.isStr(split[0]))
        self.assertTrue(isinstance(int(split[1]), int))
        self.assertTrue(util.isStr(split[2]))

    test(dm.DVH_MUTATION_FILE, cn.SPECIES_MIX_DVH)
    test(dm.MMP_MUTATION_FILE, cn.SPECIES_MIX_MMP)


if __name__ == '__main__':
    unittest.main()
