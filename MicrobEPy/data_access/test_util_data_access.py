"Utilities for data access"

import __init__
import constants as cn
import util_data_access as util_data_access


################## TEST FUNCTIONS #################
class TestFunctions(unittest.TestCase):

  def setUp(self):
    self.df = DATAFRAME.copy(deep=True)

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
      df = util_data_access.readDataModelCSV(filename)
      self.assertEqual(set(df.columns), expecteds)
    #
    test(cn.TABLE_ISOLATE)
    test("%s.csv" % cn.TABLE_ISOLATE)
    test(cn.TABLE_SCHEMAS.getSchema(cn.TABLE_ISOLATE))

  def testReadDataModelCSVCulture(self):
    if IGNORE_TEST:
      return
    df = util.readDataModelCSV(cn.TABLE_CULTURE)
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
