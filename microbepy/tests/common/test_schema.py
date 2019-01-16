from microbepy.common.schema  \
    import Schemas, ColumnSchemas, TableSchemas, FunctionalDependency
from microbepy.common import constants as cn

import copy
import collections
import unittest


COL_A = 'col_a'
COL_B = 'col_b'
COL_C = 'col_c'
DATA = range(5)
COLUMN_SCHEMAS = {COL_A: DATA, COL_B: DATA}
TABLE_A = 'table_a'
TABLE_B = 'table_b'



########################################
class TestSchema(unittest.TestCase):

  def setUp(self):
    self.schemas = Schemas()

  def testConstructor(self):
    self.assertTrue(isinstance(self.schemas.schemas, dict))

  def testGet(self):
    self.schemas.schemas = COLUMN_SCHEMAS
    self.assertEqual(self.schemas.getSchema(COL_A), DATA)

  def testValidate(self):
    self.schemas.schemas = COLUMN_SCHEMAS
    # The following should not generate an exception
    self.schemas.validate(COL_A, is_present=True)
    self.schemas.validate(COL_C, is_present=False)
    self.schemas.validate([COL_A, COL_B], is_present=True)
    #
    with self.assertRaises(ValueError):
      self.schemas.validate(COL_C, is_present=True)
    with self.assertRaises(ValueError):
      self.schemas.validate(COL_B, is_present=False)


class TestColumnSchemas(unittest.TestCase):

  def setUp(self):
    self.schemas = ColumnSchemas()

  def testConstructor(self):
    self.assertTrue(isinstance(self.schemas.schemas, dict))

  def testAdd(self):
    self.schemas.addSchema(COL_A, data_type=int)
    self.assertEqual(self.schemas.getType(COL_A), int)
    self.schemas.addSchema(COL_B)
    self.assertEqual(self.schemas.getType(COL_B), str)
    with self.assertRaises(ValueError):
      self.schemas.addSchema(COL_B)


class TestTableSchemas(unittest.TestCase):

  def setUp(self):
    self.schemas = TableSchemas()

  def testConstructor(self):
    self.assertTrue(isinstance(self.schemas.schemas, dict))

  def testAddSchema(self):
    self.schemas.column_schemas.addSchema(
        [COL_A, COL_B], data_type=int)
    self.schemas.addSchema(TABLE_A, [COL_A, COL_B], COL_A)
    self.assertEqual(self.schemas.getSchema(TABLE_A).name, TABLE_A)

  def testAddFD(self):
    self.schemas.column_schemas.addSchema(
        [COL_A, COL_B], data_type=int)
    self.schemas.addFD(COL_A, COL_B)
    self.assertEqual(self.schemas.functional_dependencies[0],
        FunctionalDependency(ind=COL_A, dep=COL_B))

  def testGetColumns2(self):
    columns = cn.TABLE_SCHEMAS.getColumns([
        cn.TABLE_MUTATION,
        cn.TABLE_ISOLATE_MUTATION_LINK,
        cn.TABLE_ISOLATE,
        cn.TABLE_GENE_DESCRIPTION,
        ])
    self.assertTrue(cn.KEY_ISOLATE in columns)


if __name__ == '__main__':
    unittest.main()
