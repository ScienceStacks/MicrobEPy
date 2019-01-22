"""
Utilities
"""

from microbepy.common import combination_iterator
from microbepy.common import config
from microbepy.common import constants as cn
from microbepy.common.equivalence_class import EquivalenceClass
from microbepy.common import schema

import collections
import math
import matplotlib.cm as cm
import numpy as np
import os
import pandas as pd
from sklearn import linear_model
import sqlite3 as sql
import sys
import types

DATA_DIR = "Data"
ALT_DATA_DIR = 'data_base'  # Directory if project is installed
REFERENCE_DIR = "reference"
DATA_MODEL_DIR = "data_model"
DATA_DIRECTORIES = set([DATA_DIR, ALT_DATA_DIR])
GIT_DIR = ".git"


PYTHON_SUBDIRECTORIES = [
    "statistics", "model", "correlation",
    "data", "plot", "search", "common",
    ]
SPECIES = [cn.SPECIES_MIX_DVH, cn.SPECIES_MIX_MMP]
TOLERANCE = 0.001
GENE_PREFIXES = ['MMP', 'DVU']
SMALL_DETERMINANT = 1e-6
TABLES_GENOTYPE_PHENOTYPE = [
      cn.TABLE_CULTURE, cn.TABLE_MUTATION, cn.TABLE_ISOLATE,
      cn.TABLE_CULTURE_ISOLATE_LINK, cn.TABLE_ISOLATE_MUTATION_LINK,
      ]

Venn = collections.namedtuple('Venn', ['both', 'only1', 'only2'])

def isNumber(obj):
  try:
    _ = float(obj)
    return True
  except:
    return False

def isStr(obj):
  """
  :return bool: True if string or unicode
  """
  type_string = str(type(obj))
  if 'str' in type_string or 'unicode' in type_string:
    return True
  else:
    return False

def toNumber(obj):
  """
  Attempts to convert object to a number.
  First an int, then a float.
  If all fail, returns np.nan.
  """
  for typ in [int, float]:
    try:
      val = typ(obj)
      return val
    except:
      pass
  return np.nan
    

def isStr(obj):
  """
  :return bool: True if string or bytes
  """
  type_string = str(type(obj))
  if 'str' in type_string or 'bytes' in type_string:
    return True
  else:
    return False

def isFloatsEqual(fp1, fp2, tolerance=TOLERANCE):
  """
  Tests for approximate equality between two floating points.
  :param float fp1:
  :param float fp2:
  :param float tolerance: a positive number less than 1
  :return bool:
  """
  if fp2 == 0:
    if fp1 == 0:
      return True
    else:
      return False
  return (fp1 - fp2)/fp2 < tolerance

def getPath(path_initial, path_last):
  """
  Constructs a path from the initial segment and the
  termination (path_list).
  :param list-of-str path_initial:
  :param str/None pat_list:
  :return str:
  """
  def makeItem(item):
    """
    Constructs a list containing either just the single item
    (if it's non-None) or an empty list.
    :param object item:
    :return list:
    """
    if item is None:
      return []
    else:
      return [item]

  path_elements = list(path_initial)
  addendum = makeItem(path_last)
  path_elements.extend(addendum)
  #
  path = path_elements[0]
  if len(path_elements) > 1:
    for ele in path_elements[1:]:
      path = os.path.join(path, ele)
  return path

def getDataModelPath(filename):
  """
  :return str path: path to data model.
  :param str filename or None: name of file in sequence data
     if None, returns the path to the directory
  """
  if cn.IS_TEST:
    return cn.TEST_PATH
  else:
    try:
      path = getIdentifiedDirectory(key_directory=GIT_DIR)
      result = getPath([path, DATA_DIR, DATA_MODEL_DIR],
        filename)
    except ValueError:
      path = getIdentifiedDirectory(key_directory=ALT_DATA_DIR)
      result = getPath([path, ALT_DATA_DIR], filename)
  return result

def getReferenceDataPath(filename):
  """
  :return str path: path to processed data
  :param str filename or None: name of file in sequence data
     if None, returns the path to the directory
  """
  root_directory = getIdentifiedDirectory()
  return getPath([root_directory, DATA_DIR, REFERENCE_DIR],
      filename)

def getSequenceDataPath(filename):
  """
  Provides a path to the sequence data.
  :param str filename or None: name of file
     if None, returns the path to the directory
  """
  return getPath([getRootDataDirectory(),
      "sequence_data"], filename)

def getRateYieldDataPath(filename):
  """
  Provides a path to the rate-yield data.
  :param str filename or None: name of file
     if None, returns the path to the directory
  """
  return getPath([getRootDataDirectory(),
      "growth_data", "rate_yield"], filename)

def getODTimeseriesDataPath(filename):
  """
  Provides a path to the rate-yield data.
  :param str filename or None: name of file
     if None, returns the path to the directory
  """
  return getPath([getRootDataDirectory(),
      "growth_data", "OD_timeseries"],  filename)

def getGeneNamesFromList(columns):
  """
  :param list-of-str columns: list of names, some of which are genes
  :return list-of-str:
  """
  names = []
  for pfx in GENE_PREFIXES:
    pfx_len = len(pfx)
    for ele in columns:
      if not isNull(ele):
        if ele[0:pfx_len] == pfx:
          names.append(ele)
  return names

def isStr(v):
  if isinstance(v, str):
    return True
  if isinstance(v, bytes):
    return True

def messageConsole(msg):
  print("*** %s" % msg)

def isNan(v):
  if isinstance(v, float):
    return np.isnan(v)
  else:
    return False

def isNull(v):
  if isNan(v):
    return True
  singleton_types = [str, int, bytes, type(None)]
  if any([isinstance(v, t) for t in singleton_types]):
    if v is None:
      return True
    if (v == "None") or (v == "none") or (v == "nan"):
      return True
  return False

def isNanInDataFrame(df, nan_columns=None):
  """
  :param pd.DataFrame df:
  :return bool, list-of-str: list of nan columns
  """
  if nan_columns is None:
    nan_columns = []
  columns = []
  for col in df.columns.tolist():
    if col in nan_columns:
      continue
    if any([isNan(x) for x in df[col].tolist()]):
      columns.append(col)
  if len(columns) == 0:
    return False, columns
  else:
    return True, columns


def replaceNan(df, columns=None, value=0):
  """
  Replaces nan values in specified columns of a dataframe.
  :param pd.DataFramed df:
  :param list-of-str columns: columns to be transformed
  :param object value:
  """
  if columns is None:
    columns = df.columns.tolist()
  for column in columns:
    new_values  = [value if isNan(x) else x for x in df[column]]
    df[column] = new_values
    
def getColorMap(num_colors=20):
  """
  :return list-of-color:
  """
  # generate data
  N = num_colors
  x = np.random.random(size=N) * 100
  y = np.random.random(size=N) * 100
  radii = np.random.random(size=N) * 1.5
  
  # get a colormap from matplotlib
  colormap =cm.get_cmap("gist_rainbow") #choose any matplotlib colormap here
  
  # define maximum and minimum for cmap
  colorspan=[40,140]
  
  # create a color channel with a value between 0 and 1
  # outside the colorspan the value becomes 0 (left) and 1 (right)
  cmap_input=np.interp(np.sqrt(x*x+y*y),colorspan,[0,1],left=0,right=1)
  
  # use colormap to generate rgb-values
  # second value is alfa (not used)
  # third parameter gives int if True, otherwise float
  A_color=colormap(cmap_input,1,True)
  
  # convert to hex to fit to bokeh
  bokeh_colors = ["#%02x%02x%02x" % (r, g, b) for r, g, b in A_color[:,0:3]]
  return bokeh_colors

def getColumnsFromDataFrame(df, columns):
  """
  :return pd.DataFrame: has the specified columns
  :raises ValueError: column not present
  """
  trues = [c in df.columns.tolist() for c in columns]
  if not all(trues):
    raise ValueError("Column not present in DataFrame")
  return pd.DataFrame(df[columns])

def changeColumnValues(df, func):
  """
  Change values in the columns.
  :param DataFrame df:
  :param Function func:
     Inputs: column (str), value (float)
     Output: float
  """
  for col in df.columns:
    values = df[col].tolist()
    new_values = [func(col, v) for v in values]
    df[col] = new_values

def setNoneList(values):
  """
  :param list/None values:
  :return list
  """
  if values is None:
    return []
  else:
    return values

def getColumnType(column):
  """
  Finds the type for a column
  :param str column:
  :return type:
  """
  for column_schema in cn.TABLE_SCHEMAS.column_schemas.getSchemas():
    if column == column_schema.name:
      return column_schema.data_type
  raise RuntimeError("Column %s is not typed" % column)

def cleanDF(df, is_reset_index=False):
  """
  Removes duplicates and stray columns.
  Removes all columns of the form "index_".
  :param pd.DataFrame df:
  """
  df.drop_duplicates(inplace=True)
  df = pruneNullRows(df)
  if is_reset_index:
    resetIndex(df)
  remove_column = "%s_" % cn.INDEX
  for column in df.columns:
    if (column.find(remove_column) == 0)  \
       or (column == cn.INDEX):
      del df[column]
  return df

def resetIndex(df):
  """
  Resets the index handling stray columns.
  :param pd.DataFrame df:
  """
  LEVEL_0 = 'level_0'
  columns = df.columns
  if LEVEL_0 in columns:
    del df[LEVEL_0]
  df.reset_index(inplace=True)
  if LEVEL_0 in columns:
    del df[LEVEL_0]

def standardize(df, columns=None):
  """
  Standardizes a numeric column in a dataframe.
  :param pd.DataFrame df:
  :param list-of-str columns
  """
  if columns is None:
    columns = df.columns
  for column in columns:
    std = np.std(df[column])
    avg = np.mean(df[column])
    if np.isclose(std, 0):
      df[column] = 0
    else:
      df[column] = (df[column] - avg)/std

def deleteColumns(df, columns, is_drop_duplicates=True):
  """
  Deletes columns from the DataFrame, if they are present.
  :param pd.DataFrame df:
  :param list-of-str columns:
  :param bool is_drop_duplicates:
  """
  for column in columns:
    if column in df.columns:
      del df[column]
  if is_drop_duplicates:
    df.drop_duplicates(inplace=True)

def trimDF(df, keep_columns=None, delete_columns=None):
  """
  Keeps/deletes columns in a dataframe.
  Exactly one of keep_columns/delete_columns can be non-None.
  :param pd.DataFrame df:
  :param list-of-str keep_columns: Columns to keep
  :param list-of-str delete_columns: Columns to delete
  :return pd.DataFrame:
  """
  if keep_columns is None and delete_columns is None:
    raise ValueError("Invalid parameters.")
  if not keep_columns is None and not delete_columns is None:
    raise ValueError("Invalid parameters.")
  if keep_columns is not None:
    df_result = df[keep_columns].copy(deep=True)
  if delete_columns is not None:
    df_result = df.copy(deep=True)
    deleteColumns(df_result, delete_columns, is_drop_duplicates=True)
  df_result.drop_duplicates(inplace=True)
  return df_result

def makeNullRow(df, null=np.nan):
  """
  Creates a row of null values to the dataframe.
  :param object null: value in row
  :return dict: row
  """
  row = {}
  for column in df.columns.tolist():
    row[column] = null
  return row

def addNullRow(df, null=np.nan):
  """
  Adds a row of null values to the dataframe.
  :param object null: value in row
  :return pd.DataFrame:
  """
  return df.append(makeNullRow(df), ignore_index=True)

def getDBPath():
  """
  There are 3 possibilities for the database path
  1. in the Data directory at the root of the microbepy project
  2. in the Data directory at the root of a containing project
  3. it is pointed by the cn.CONFIG_FILE
  """
  yaml_dict = config.setup(create_config=False)
  if yaml_dict[cn.SQLDB_PATH_NAME] is not None:
    path = yaml_dict[cn.SQLDB_PATH_NAME]
  else:
    filename = cn.SQLDB_FILE
    try:
      srcdir = getIdentifiedDirectory(key_directory=ALT_DATA_DIR)
    except:
      srcdir = None
    if srcdir is not None:
      path = os.path.join(srcdir, ALT_DATA_DIR)
      path = os.path.join(path, filename)
    else:
      path = getDataModelPath(filename)
  if path is None:
    raise ValueError("***You must setup %s" % cn.CONFIG_FILE)
  return path

def getDBConnection():
  path = getDBPath()
  try:
    conn = sql.connect(path)
  except:
    raise ValueError("Invalid path to DB: %s" % path)
  return conn

def getDuplicates(values):
  """
  :param list-of-object values:
  :return list-of-object:
  """
  result = list(values)
  [result.remove(x) for x in set(values)]
  return result

def typeDF(df, columns=None):
  """
  Ensures that column values are of the correct type.
  :param list-of-str columns:
  :return pd.DataFrame:
  """
  if columns is None:
    columns = list(df.columns)
  for column in columns:
    typ = getColumnType(column)
    values = df[column]
    df[column] = [np.nan if isNull(v) else typ(v) for v in values]
  return df

def removeDuplicateColumns(df):
  """
  Removes columns that have a duplicate name.
  :return pd.DataFrame:
  """
  duplicates = getDuplicates(df.columns)
  done = False
  idx = 0
  df_result = df.copy()
  additions_dict = {}
  while not done:
    if idx >= len(df_result.columns):
      done = True
      break
    column = df_result.columns[idx]
    if column in duplicates:
      df1 = df_result[column]
      values = df1.iloc[:,1]
      del df_result[column]
      duplicates.remove(column)
      additions_dict[column] = values
    else:
      idx += 1
  df_add = pd.DataFrame(additions_dict)
  df_result = pd.concat([df_result, df_add], axis=1, sort=True)
  return df_result

def makeVenn(group1, group2):
  """
  Constructs the Venn regions for two groups:
    both - overlap between the groups
    only1 - only n group 1
    only2 - only n group 2
  :param list-of-object group1:
  :param list-of-object group2:
  :return Venn:
  """
  set1 = set(group1)
  set2 = set(group2)
  return Venn(
      both=set1.intersect(set2),
      only1=set1.difference(set2),
      only2=set2.difference(set1),
      )

def selNonNull(values1, values2):
  """
  Selects the non-null value of the two lists in its position.
  :param list-of-object values1:
  :param list-of-object values2:
  :param list-of-object result:
  """
  if len(values1) != len(values2):
    raise ValueError("Inputs must have the same length.")
  pairs = zip(values1, values2)
  result = []
  for x,y in pairs:
    if x == y:
      result.append(x)
    elif isNull(x):
      result.append(y)
    elif isNull(y):
      result.append(x)
    else:
      msg = "Inputs must be identical if in the same position and non-null."
      raise ValueError(msg)
  return result
  
def mergeRowsColumns(df1, df2, merge_column):
  """
  Merges the two dataframes, combining both rows and columns.
  Column values are appended where the dataframes have the same columns.
  nan values are added to columns where the column is not present
  in a dataframe.
  :param pd.DataFrame df1:
  :param pd.DataFrame df2:
  :param str merge_column: column on which the merge is done
  :return pd.DataFrame: Has columns df1.columns + df2.columns
  """
  LEFT = "__left"
  RIGHT = "__right"
  df_1 = removeDuplicateColumns(df1)
  cleanDF(df_1)
  df_2 = removeDuplicateColumns(df2)
  cleanDF(df_2)
  df_result = df_1.merge(df_2, on=merge_column, how='outer',
      suffixes=(LEFT, RIGHT))
  # Merge the overlapping columns
  left_overlaps = [c for c in df_result.columns
                   if c[-len(LEFT):] == LEFT]
  left_overlaps.sort()
  right_overlaps = [c for c in df_result.columns
                   if c[-len(RIGHT):] == RIGHT]
  right_overlaps.sort()
  pairs = zip(left_overlaps, right_overlaps)
  for left, right in pairs:
    values = selNonNull(df_result[left], df_result[right])
    pos = left.find(LEFT)
    column = left[0:pos]
    del df_result[left]
    del df_result[right]
    df_result[column] = values
  # Finalize result
  cleanDF(df_result)
  return df_result

def readSQL(cmd):
  """
  Creates a dataframe for the SQL query.
    1. Duplicate column names (ending with ':n') are merged.
    2. Deletes columns that begin with 'index'
    3. None values have a consistent representation
  :param str/TableSchema cmd: SQL query command or TableSchema,
      if want the entire table
  :return pd.DataFrame:
  """
  SEP = ":"
  if isinstance(cmd, schema.TableSchema):
    sql_cmd = "SELECT * FROM %s" % cmd.name
  else:
    sql_cmd = cmd
  conn = getDBConnection()
  df_result = pd.read_sql(sql_cmd, conn)
  conn.close()
  done = False
  while not done:
    duplicates = [c for c in df_result.columns if c.count(SEP) > 0]
    if len(duplicates) == 0:
      done = True
      break
    column = duplicates[0]
    pos = column.find(SEP)
    column = column[:pos]
    columns = [c for c in df_result.columns if c.count(column) > 0]
    # Delete index columns
    if column in ['level_0', cn.INDEX]:
      for col in columns:
        del df_result[col]
    # Merge other columns
    else:
      col_first = columns[0]
      values = df_result[col_first].tolist()
      del df_result[col_first]
      columns = columns[1:]
      for col in columns:
        try:
          values = selNonNull(values, df_result[col].tolist())
        except:
          import pdb; pdb.set_trace()
        del df_result[col]
      df_result[column] = values
  df_result = unifyNullValues(df_result)
  return df_result

def unifyNullValues(df_value):
  """
  Makes all null df_value cn.NONE.
  :param pd.DataFrame df_value:
  :return pd.DataFrame:
  """
  return df_value.fillna(value=np.nan)
  #return df_value.applymap(lambda v: cn.NONE if isNull(v) else v)

def pruneNullRows(df):
  """
  Removes rows that are all nulls.
  :param pd.DataFrame df:
  This is done in place to avoid storage problems with large dataframes.
  :return pd.DataFrame:
  """
  return df.dropna(axis=0, how='all')

def pruneRowsWithNullColumns(df, columns):
  """
  Deletes rows where there is a null value in any of a list
  of columns.
  :param pd.DataFrame df:
  :param list-of-str colums:
  :return pd.DataFrame:
  """
  def check(row):
    return not any([isNull(row[c]) for c in columns])
  #
  sel = df.apply(check, axis=1)
  return pd.DataFrame(df.loc[sel])

def coerceDF(df):
  """
  Coerces the columns to their type if the type is known.
  :param pd.DataFrame df:
  :return pd.DataFrame:
  """
  df_result = df.copy(deep=True)
  for column in df_result.columns:
    try:
      schema = cn.TABLE_SCHEMAS.column_schemas.getSchema(column)
      if schema.data_type in [float, int, bool]:
        df_result[column] = pd.to_numeric(df_result[column])
    # Get an exception if the column type is unknown
    except ValueError:
      pass
  return df_result

def appendUnique(values, value):
  """
  Adds value to values if it is not already present.
  :param list-of-object values:
  :param object value:
  :return list-of-object:
  """
  new_values = list(values)
  if not value in new_values:
    new_values.append(value)
  return new_values

def mergeLeft(df1, df2, merge_column):
  """
  Does a left join properly handling null values (which pandas does not do).
  :param pd.DataFrame df1: left dataframe
  :param pd.DataFrame df2:
  :param str merge_column: column on which the merge is done
  :return pd.DataFrame:
  Assumes that the only overlap in column names is the merge_column
  """
  overlaps = set(df1.columns).intersection(df2.columns)
  if overlaps != set([merge_column]):
    raise ValueError ("Must only have the merge column in common!")
  # Find the null values on the left
  sel = [not x for x in df1[merge_column].isnull()]
  df_left = df1[sel]
  # Merge non-null values
  df_result = df_left.merge(df2, on=merge_column, how='left')
  # Include the omitted rows as null columns
  sel = [x for x in df1[merge_column].isnull()]
  df_null = df1[sel].copy()
  for column in df2.columns:
    if column != merge_column:
      df_null[column] = cn.NONE
  # Add the missing rows
  df_result = pd.concat([df_result, df_null], sort=True)
  # Finalize result
  cleanDF(df_result)
  return df_result

def unpivot(df, label_column='label', value_column='value'):
  """
  Converts a pivoted dataframe into one with a column
  of labels.
  :param pd.DataFrame: DF with values in cells and value
                       labels in columns
  :param str label_column: name of the label column on output
  :param str value_column: name of the value column on output
  :return pd.DataFrame: label_column, value_column
  """
  return pd.melt(df, var_name=label_column, value_vars=df.columns,
      value_name=value_column)

def makeMatrix(df, 
      row_name=cn.KEY_ISOLATE,
      column_name=cn.KEY_MUTATION, 
      value_name=cn.COUNT,
      default_value=0.0,
      ):
  """
  Creates a data matrix that consists of only a column of
  row names, a column of column name, and a column of values.
  :param pd.DataFrame df:
  :param str row_name: name of the column in df used for row values
  :param str column_name: name of column whose values are used
                          as column names in the matrix
  :param str value_name: name of the column from which values
                         are obtained.
  :param float default_value: used for missing values
  """
  df_sub = df[[row_name, column_name, value_name]].copy()
  df_sub.drop_duplicates(inplace=True)
  sel = df_sub.apply(
     lambda r: (not isNull(r[column_name]))
     and (not isNull(r[row_name])),
     axis=1
     )
  df_sub = df_sub.loc[sel]
  df_result = df_sub.pivot_table(index=row_name, columns=column_name, values=value_name)
  df_result = df_result.applymap(lambda x: 0 if isNull(x) else x)
  return df_result

def set2list(aSet):
  return [x for x in aSet]

def removeDuplicatesFromList(values):
  """
  Trims the list to remove duplicates. List elements
  may be a string or a list of strings or int.
  Order is not preserved.
  :param list-of-str or list-of-list-of-other values:
  :return list-of-inputtype:
  """
  if isStr(values[0]) or isNumber(values[0]):
    new_values = set(values)
    return set2list(new_values)
  else:
    df = pd.DataFrame(values)
    df = df.drop_duplicates()
    results = []
    for _,row in df.iterrows():
      results.append([v for v in row.to_dict().values()])
    return results

def getParentClass(cls):
  return cls.__bases__[0]

def normalizedMean(values, max_std):
  """
  Computes the mean of values adjusted by their standard deviation.
  :param list-of-number values:
  :param float max_std:
  :return float:
  """
  trimmed_values = trimExtremeValues(values, max_std)
  std = np.std(trimmed_values)
  if np.isclose(std, 0):
    return 0.0
  return np.mean([v/std for v in trimmed_values])

def trimExtremeValues(values, max_std):
  """
  Returns a list that has values no large than max_std.
  :param list-of-number values:
  :param float max_std:
  :return list-of-number:
  """
  std = np.std(values)
  if np.isclose(std, 0):
    return values
  normalized_values = [v/std for v in values]
  pairs = zip(values, normalized_values)
  return [v for v,z in pairs if abs(z) <= max_std]

def selectRows(df, constraints):
  """
  Finds a subset of the original DataFrame.
  :param pd.DataFrame df:
  :param list-of-BooleanFunction constraints: 
      A boolean function of a row of df_culture_isolate_mutation.
  :return pd.DataFrame:
  """
  constraints = setNoneList(constraints)
  def checkConstraint(row):
    is_included = True
    for func in constraints:
      if not func(row):
        is_included = False
        break
    return is_included
  #
  sel = df.apply(checkConstraint, axis=1)
  df_result = df.loc[sel].copy(deep=True)
  df_result.drop_duplicates(inplace=True)
  return df_result

def isColinear(df_X):
  """
  Determines if there is a colinearity in the X matrix.
  :param pd.DataFrame df_X: X matrix in linear regression
  :return bool:
  """
  pseudo_inverse = np.dot(df_X.transpose(), df_X)
  det = np.linalg.det(pseudo_inverse)
  is_colinear = np.abs(det) < SMALL_DETERMINANT
  return is_colinear

def findRowColumn(df, predicate):
  """
  Returns index, column tuples whose value satisfies a predicate.
  :param pd.DataFrame df:
  :param UnaryBooleanFunction predicate:
  :return list-of-tuple:
  """
  result = []
  for index in df.index:
    for column in df.columns:
      if predicate(df.loc[index, column]):
        result.append((index, column))
  return result

def aggregateColumns(df, columns, aggregateFunc, sep="--"):
  """
  Creates a new column that that is an aggregation of existing columns.
  The old column is deleted. The new column is named
  as the concatenation of old columns separated by the separator.
  :param pd.DataFrame df:
  :param list-str columns:
  :param Function aggregateFunc: arg is dataframe; returns Series
  :return str: name of the new column
  """
  df_sub = df[list(columns)]
  merged = aggregateFunc(df_sub)
  for col in columns:
    del df[col]
  str_columns = [str(c) for c in columns]
  new_column = sep.join(str_columns)
  df[new_column] = merged
  return new_column

def findCorrelatedColumns(df, min_corr=0.99):
  """
  Forms groups of columns based on their correlations.
  :param pd.DataFrame df: Dataframe whose columns are analyzed
  :param float min_corr: minimum correlation for two mutations to be equivalent
  :return list-list-columns: List of columns that are correlated.
  """
  df_corr = df.corr()
  relation = lambda m1, m2: df_corr.loc[m1, m2] >= min_corr
  equiv = EquivalenceClass(df, relation)
  equiv.do()
  equiv.validate()  # Verify that the relations hold for groups
  classes = [c for c in equiv.classes if len(c) > 1]
  return classes

def isIterable(values):
  return "__iter__" in dir(values)

def selectColumns(df, constraint):
  """
  Selects a subset of the columns in a dataframe.
  :param pd.DataFrame df:
  :param BooleanFunction constraint: Boolean Function of column name
  :return pd.DataFrame:
  """
  columns = [c for c in df.columns if constraint(c)]
  return df[columns].copy()

def correlateWithPredictors(df_X, df_y):
  """
  Constructs the correlation of each predictor with a dependent
  variable.
  :param pd.DataFrame df_X: columns of predictor variables
  :param pd.DataFrame df_y: dependent variable
  :return pd.DataFrame: cn.VALUE, indexed by df_column, ordered
                        by descending absolute value of correlation
  """
  df = df_X.copy()
  df[cn.DEPVAR] = df_y.iloc[:,0]
  df_corr = df.corr()
  df_corr = df_corr.applymap(lambda v: abs(v))
  df_result = pd.DataFrame(df_corr[cn.DEPVAR])
  df_result = df_result.drop([cn.DEPVAR])
  df_result = df_result.sort_values(cn.DEPVAR, ascending=False)
  return df_result

def appendWithColumnUnion(df1, df2, fill_value=0.0):
  """
  Appends the two dataframes, adding columns if needed.
  :param pd.DataFrame df1:
  :param pd.DataFrame df2:
  :param float default_value: fill value used for added columns
  """
  def addColumns(df, columns):
    for col in columns:
      df[col] = fill_value
  #
  columns1 = set(df1.columns)
  columns2 = set(df2.columns)
  add_columns1 = columns2.difference(columns1)
  add_columns2 = columns1.difference(columns2)
  #
  df1_new = df1.copy()
  df2_new = df2.copy()
  #
  addColumns(df1_new, add_columns1)
  addColumns(df2_new, add_columns2)
  df_result = df1_new.append(df2_new)
  #
  df_result.index = range(len(df_result))
  return df_result

def makeProjectedPredictor(df_X, df_x):
  """
  Computes the residuals of df_x regressed on df_X and
  returns a matrix that augments df_X with these residuals.
  :param pd.DataFrame df_X:
  :param pd.DataFrame df_x:
  :return pd.DataFrame:
  """
  MIN_SCORE = 0.05  # Minim score to use the residuals
  model = linear_model.LinearRegression(fit_intercept=False, 
      copy_X=True)
  if isinstance(df_x, pd.Series):
    df_x = pd.DataFrame(df_x)
  name = df_x.columns[0]
  if isColinear(df_X):
    df_res = df_x
  else:
    model.fit(df_X, df_x)
    df_res = pd.DataFrame({
        name: [v[0] for v in model.predict(df_X)],
        })
    score = model.score(df_X, df_x)
    if score < MIN_SCORE:
      df_res = df_x
    else:
      df_res = pd.DataFrame(df_x[name] - df_res[name])
      df_res = df_res.applymap(lambda v: 0.0 if np.isnan(v) else v)
  df_result = df_X.copy()
  df_result[name] = df_res[name]
  return df_result

def makeProjectedPredictors(df_X):
  """
  Computes the residuals of the regression of successive columns
  on preceeding columns.
  :param pd.DataFrame df_X:
  :return pd.DataFrame:
  """
  columns = df_X.columns.tolist()
  df_result = pd.DataFrame(df_X[columns[0]])
  columns = columns[1:]
  for col in columns:
    df = makeProjectedPredictor(df_result, df_X[col])
    df_result[col] = df[col]
  return df_result

def extendBin(bnum, length):
  """
  :param str bnum: binary number or integer
  :param int length:
  :return list-str: binary number with leading zeros
  """
  if not isinstance(bnum, str):
    bnum = bin(bnum)
  result = list(bnum)[2:]
  while len(result) < length:
    result.insert(0, '0')
  return result

def addBinColumn(df):
  """
  For a dataframe with columns that are binary values,
  adds a column that is the binary number for that column.
  On return, df has a column cn.VALUE that is the decimal
  number for the binary value.
  :param pd.DataFrame df:
  """
  binaries = []
  columns = df.columns
  dff = df.copy()
  for col in dff.columns:
    dff = dff.rename(columns={col: str(col)})
  for _, row in dff.iterrows():
    binary = '0b'
    for col in dff.columns:
      binary = binary + str(int(row[col]))
    binaries.append(int(binary, 2))
  df[cn.VALUE] = binaries

def getFirstColumn(df):
  """
  :param pd.DataFrame df:
  :return list-object: first column of df
  """
  return df[df.columns[0]].tolist()

def nCr(n, r):
  """
  Number of combinations of a set of size n of distinct objects
  choosing r objects.
  :param int n:
  :param int r:
  :return int:
  """
  return math.factorial(n) / math.factorial(n-r) / math.factorial(r)

def findKey(df, required_columns=None, max_size=4):
  """
  Finds a key for the dataframe.
  :param pd.DataFrame df:
  :param list-str required_columns: Columns to include in the key
  """
  length = len(df.drop_duplicates())
  required_columns = setNoneList(required_columns)
  columns = df.columns.tolist()
  for col in required_columns:
    columns.remove(col)
  combinator = combination_iterator.CombinationIterator(columns,
      max_size)
  keys = []
  for new_columns in combinator:
    new_columns.extend(required_columns)
    if len(df[new_columns].drop_duplicates()) == length:
      keys = new_columns
      break
  return keys
    
def makeGeneDescriptionDF(ggene_ids):
  """
  Creates a dataframe of descriptions.
  :param list-of-str ggene_ids: values in cn.GGENE_ID
  :return pd.DataFrame: cn.GGENE_ID, GENE_DESC
  """
  quoted_list = str(["%s" % str(g) for g in ggene_ids])
  in_list = quoted_list.replace("[", "")
  in_list = in_list.replace("]", "")
  query = '''
      select distinct gene_id, ggene_id, gene_desc
      from genotype 
      where gene_id = ggene_id and ggene_id in (%s)
      ''' % in_list
  df = readSQL(query)
  return df

def getIsolatesFromIndices(indices):
  """
  Extracts the isolates from the indices of a df_X.
  :param pandas.index indices: 
      cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP
  :return dict: keyed by cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP
    values correspond to rows element in the index
  """
  keys = [n for n in indices.names]
  result = {}
  for idx, key in enumerate(keys):
    result[key] = [v[idx] for v in indices.values]
  return result

def makeDataframeFromXlsx(path, sheet_no=0):
  """
  Creates a dataframe from an xlsx file in a data directory
  :param str path: Path to data files
  :param int sheet_no: Number of the worksheet
  """
  data = pd.ExcelFile(path)
  return data.parse(sheet_no)

def getRootDataDirectory():
  return getPath([getIdentifiedDirectory(), DATA_DIR], None)

def getIdentifiedDirectory(key_directory=None):
  """
  The root directory is the root of the enclosing project
  (since microbepy is intended to be a submodule).
  :return str: path to top folder of enclosed project
  """
  if key_directory is None:
    directories = DATA_DIRECTORIES
  else:
    directories = set([key_directory])
  curdir = os.getcwd()
  paths = []
  while len(curdir) > 1:
    paths.append(curdir)
    curdir = os.path.split(curdir)[0]
  paths.reverse()
  for path in paths:
    if len(directories.intersection(os.listdir(path))) > 0:
      return path
  raise ValueError("%s not found." % key_directory)

def makeNormalizedData(df_denormalized):
  """
  Creates a set of normalized dataframes that correspond to
  the denomralized data for genotype_phenotype data.
  :param pd.DataFrame df_denormalized: 
      See constants.TABLE_CULTURE_ISOLATE_MUTATION
  :return NormalizedData:
  """
  normalized_data = NormalizedData()
  for table_name in TABLES_GENOTYPE_PHENOTYPE:
    schema = cn.TABLE_SCHEMAS.getSchema(table_name)
    df = df_denormalized[schema.columns].copy()
    normalized_data[schema.name] = df.drop_duplicates()
  return normalized_data

def makeDenormalizedDF(normalized_data):
  """
  Creates a denormalized dataframe by joining the normalized tables.
  :param NormalizedData normalized_data:
  :return pd.DataFrame: has schema of TABLE_CULTURE_ISOLATE_MUTATION
  """
  def mergeDF(table1, table2, table_link, columns):
    df1 = normalized_data[table1].merge(normalized_data[table_link],
        on=columns[0], how='inner')
    df2 = df1.merge(normalized_data[table2],
        on=columns[1], how='inner')
    df = df2.drop_duplicates()
    return df
  #
  df_im = mergeDF(cn.TABLE_MUTATION, cn.TABLE_ISOLATE,
      cn.TABLE_ISOLATE_MUTATION_LINK,
      [cn.KEY_MUTATION, cn.KEY_ISOLATE])
  df_ci = normalized_data[cn.TABLE_CULTURE].merge(
      normalized_data[cn.TABLE_CULTURE_ISOLATE_LINK], 
      on=cn.KEY_CULTURE, how="inner")
  df = df_im.merge(df_ci, on=cn.KEY_ISOLATE, how="inner")
  return df
    


################# CLASSES ####################
class NormalizedData(dict):
  """
  Container of tables for Noramlized Data.
  Has keys for TABLE_PHENOTYPE_GENOTYPE
  """
  pass
