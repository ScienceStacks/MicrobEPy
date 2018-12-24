"Utilities for data access"

import __init__
import constants as cn
import schema
import util
import util_data_access as util_data_access

import unittest

IGNORE_TEST = False


################## TEST FUNCTIONS #################
class TestFunctions(unittest.TestCase):

  def testMakeCSVFilename(self):
    if IGNORE_TEST:
      return
    expected = "isolate.csv"
    self.assertEqual(util_data_access.makeCSVFilename(cn.TABLE_ISOLATE),
        expected)
    self.assertEqual(
        util_data_access.makeCSVFilename(
        cn.TABLE_SCHEMAS.getSchema(cn.TABLE_ISOLATE)),
        expected)

  def testReadDataModelCSV(self):
    def test(filename):
      expecteds = set(
          cn.TABLE_SCHEMAS.getSchema(cn.TABLE_ISOLATE).columns)
      try:
        df = util_data_access.readDataModelCSV(filename)
      except FileNotFoundError:
        if isinstance(filename, schema.TableSchema):
          filename = filename.name
        print("CSV file is not present for %s" % filename)
        return
      self.assertEqual(set(df.columns), expecteds)
    #
    test(cn.TABLE_ISOLATE)
    test("%s.csv" % cn.TABLE_ISOLATE)
    test(cn.TABLE_SCHEMAS.getSchema(cn.TABLE_ISOLATE))

  def testReadDataModelCSVCulture(self):
    if IGNORE_TEST:
      return
    try:
      df = util_data_access.readDataModelCSV(cn.TABLE_CULTURE)
    except FileNotFoundError:
      print("CSV file is not present")
      return
    nulls = [r for _,r in df.iterrows() 
             if util.isNull(r[cn.KEY_CULTURE])]
    self.assertTrue(len(nulls) == 0)
 
  def testMakeDataframeFromXlsx(self):
    if IGNORE_TEST:
      return
    test_file = util.getSequenceDataPath("dvh_mutations.xlsx")
    df = util.makeDataframeFromXlsx(test_file)
    

if __name__ == '__main__':
    unittest.main()
