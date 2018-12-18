"""
Represents the Schemas of Tables and Columns.
Each class creates a singleton instance.
This module cannot depend on any other module in coevolution to avoid
circular references.
"""

import __init__
from collections import namedtuple
  
ColumnSchema = namedtuple('ColumnSchema', ['name', 'data_type'])
TableSchema = namedtuple('TableSchema', [
    'name',     # Table name
    'columns',  # columns in the table
    'key',       # join column
    ])
FunctionalDependency = namedtuple('FunctionalDependency',
  ['ind', 'dep'])

#################################
def assignList(value):
  if isinstance(value, list):
    result = value
  else:
    result = [value]
  return result


#################################
class Schemas(object):
  """ Abstract class for common properties and methods."""

  def __init__(self):
    self.schemas = {}

  def getAll(self):
    return self.schemas

  def addSchema(self, value):
    """
    :param Schemas or list-of-Schemas value:
    """
    schemas = assignList(value)
    for schema in schemas:
      self.validate(schema.name, is_present=False)
      self.schemas[schema.name] = schema

  def getSchema(self, name):
    self.validate(name, is_present=True)
    return self.schemas[name]

  def getSchemaNames(self):
    return list(self.schemas.keys())

  def getSchemas(self):
    return list(self.schemas.values())

  def validate(self, value, is_present=None):
    """
    :param str or list-of-str value:
    :param bool is_present: expect to be present
    :raises ValueError: if value not in schema
    """
    if is_present is None:
      raise RuntimeError("Must specify 'is_present'")
    names = assignList(value)
    for name in names:
      if name in list(self.schemas.keys()):
        if is_present:
          return
        raise ValueError("%s is present but not should be." % name)
      else:
        if is_present:
          raise ValueError("%s not present but should be." % name)
        return


##########################################
class ColumnSchemas(Schemas):

  def addSchema(self, names, data_type=str):
    """
    :param list-of-str or str value
    :param Type data_type:
    """
    names = assignList(names)
    for name in names:
      super(self.__class__, self).addSchema(
          ColumnSchema(name=name, data_type=data_type))

  def getType(self, name):
    return self.schemas[name].data_type


##########################################
class TableSchemas(Schemas):

  def __init__(self):
    super(self.__class__, self).__init__()
    self.column_schemas = ColumnSchemas() # Column schemas
    self.functional_dependencies = []

  def addSchema(self, name, columns, key):
    """
    :param str name:
    :param list-of-str columns:
    :param list-of-str or str key:
    Columns must already exists
    """
    self.column_schemas.validate(columns, is_present=True)
    self.column_schemas.validate(key, is_present=True)
    super(self.__class__, self).addSchema(
        TableSchema(name=name, columns=columns, key=key))

  def addFD(self, ind_column, dep_column):
    """
    Specifies that value of ind_column -> value of dep_column
    :param str ind_column:
    :param str dep_column:
    """
    self.column_schemas.validate(
        [ind_column, dep_column], is_present=True)
    self.functional_dependencies.append(FunctionalDependency(
        ind=ind_column, dep=dep_column))

  def getColumns(self, table_names):
    table_names = assignList(table_names)
    result = []
    for name in table_names:
      table = self.getSchema(name)
      result.extend(table.columns)
    result = list(set(result))
    return result
      
