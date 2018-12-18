"""
Provides convenient access to columns in the data model with
pandas DataFrames.

Three levels of API are provided:

  Instance variables (lowest level). Provides a DataFrame instantion
  of the CSV files in the data model. Dataframes are in self.dfs.

  Issue a SQL command against the convolution database (cn.DB_NAME).
  See readSQL().

  Projections of the join of culture, isolate, mutation
  See  makeDF().  Creates a dataframe consisting
  of an arbitrary set of columns from the three tables.

__init__
  1. Ensure that tables read have nan values as Nulls
For makeDF:
  1. Check for columns in GENE_DESCRIPTION; do "outer" join if needed.
"""


import __init__
import constants as cn
import isolate as iso
import util

import copy
import os
import numpy as np
import pandas as pd


#############################################
# Helper Functions
############################################


#############################################
# API Class
############################################
class Api(object):
  """
  Provides access to the data model at 3 levels,
  as described above.
  """

  def __init__(self, db_name=cn.DB_NAME):
    """
    :param str db_name: name of the database
    """
    # Dataframes from CSV files
    self.dfs = {}
    self.schemas = {}
    for schema in cn.TABLE_SCHEMAS.getSchemas():
      self.schemas[schema.name] = schema
      self.dfs[schema.name] = util.readSQL(schema)

  def makeDF(self, columns,
      aux_columns=None,
      constraints=None, effects=None, 
      species_mixes=None,
      no_null_columns=[cn.KEY_MUTATION, 
          cn.KEY_ISOLATE, cn.KEY_CULTURE, cn.GENE_ID, cn.GGENE_ID]):
    """
    Provides any combination of columns in TABLE_CULTURE_ISOLATE_MUTATION.
    Selects non-duplicate rows.
    Assumes that true keys (KEY_MUATION, KEY_ISOLATE, KEY_CULTURE)
    should be non-null.
    Coerces columns to their data types, although this cannot always happen
    because of NULL values.
    :param list-of-str columns: names of the columns desired. If None, then all.
        requested_columns,
    :param list-of-str aux_columns: names of the columns needed to compute
        constraints, but these columns are not returned.
    :param list-of-BooleanFunction constraints: 
        A boolean function of a row of df_culture_isolate_mutation.
    :param list-of-str effects: a constraint for the mutation effects considered
    :param list-of-str species_mix: a constraint for the species mixes considered
    :param list-of-str no_none_columns: list of columns that should no contain a "none"
    :return pd.DataFrame:
    """
    def isNoNullColumns(row):
      for column in no_null_columns:
        if column in row.index:
          if util.isNull(row[column]):
            return False
      return True
    #
    constraints = util.setNoneList(constraints)
    aux_columns = util.setNoneList(aux_columns)
    # Adjust the columns that should not contain null values
    no_null_columns = list(set(no_null_columns).intersection(columns))
    for fd in cn.TABLE_SCHEMAS.functional_dependencies:
      if (fd.ind in no_null_columns) and (fd.dep in no_null_columns):
          no_null_columns.remove(fd.dep)
    # Construct the complete set of columns to consider
    all_columns = list(columns)
    all_columns.extend(aux_columns)
    # Adjust aux_columns to include inferred constraints
    all_columns = util.appendUnique(all_columns, cn.EFFECT)
    all_columns = util.appendUnique(all_columns, cn.SPECIES_MIX)
    # Select the columns
    df_initial = None
    cim_columns = list(set(all_columns).intersection(
        self.schemas[cn.TABLE_CULTURE_ISOLATE_MUTATION].columns))
    df_initial = self.dfs[
        cn.TABLE_CULTURE_ISOLATE_MUTATION][cim_columns].copy(deep=True)
    # Handle the row constraints
    if effects is not None:
      effects_str = [str(x) for x in effects]
      # TODO: Why didn't lambda statement work?
      def checkEffect(row):
        return str(row[cn.EFFECT]) in effects_str
      #
      #constraints.append(lambda r: str(r[cn.EFFECT]) in effects_str)
      constraints.append(checkEffect)
    if species_mixes is not None:
      species_mixes_str = [str(x) for x in species_mixes]
      constraints.append(
          lambda r: str(r[cn.SPECIES_MIX]) in species_mixes_str)
    constraints.append(isNoNullColumns)
    df_result = self.selectRows(df_initial, constraints)
    # Trim the result
    df_result = df_result[columns].copy(deep=True)
    df_result.drop_duplicates(inplace=True)
    util.cleanDF(df_result)
    # Coerce column values
    df_result = util.coerceDF(df_result)
    #
    return df_result

  def selectRows(self, df, constraints):
    """
    Convenience method
    """
    return util.selectRows(df, constraints)

  def readSQL(self, sql_cmd):
    """
    Queries the convolution database using the sql command.
    :param str sql_cmd:
    :return pd.DataFrame:
    """
    df = util.readSQL(sql_cmd)
    return util.coerceDF(df)  # Adjust values to their types
