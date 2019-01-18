from microbepy.common import helpers
from microbepy.common import util
from microbepy.common import constants as cn

import copy
import collections
import numpy as np
import pandas as pd
import random
import os
import unittest


DIR_SEPARATOR = "/"
IGNORE_TEST = False
COL_A = 'a'
COL_B = 'b'
COL_C = 'c'
COL_D = 'd'
SIZE = 10
DATA = list(range(10))
COLUMNS = [COL_A, COL_B, COL_C]
DATAFRAME = pd.DataFrame()
for column in COLUMNS:
  DATAFRAME [column] = list(range(SIZE))

################## HELPERS #################
makeData = lambda: [random.normalvariate(0, 1) for _ in range(SIZE)]


################## TEST FUNCTIONS #################
class TestFunctions(unittest.TestCase):

  def setUp(self):
    self.df = DATAFRAME.copy(deep=True)

  def testGetProjectDirectory(self):
    if IGNORE_TEST:
      return
    result = util.getIdentifiedDirectory()
    split_result = result.split(DIR_SEPARATOR)
    cur_path = os.getcwd()
    split_cur_path = cur_path.split(DIR_SEPARATOR)
    length = len(split_result)
    expected = split_cur_path[:length]
    self.assertEqual(split_result, expected)

  def testGetPaths(self):
    if IGNORE_TEST:
      return
    funcs = [util.getSequenceDataPath,
        util.getRateYieldDataPath, util.getODTimeseriesDataPath]
    for func in funcs:
      result = func("dummy")
      self.assertTrue(isinstance(result, str))
    self.assertTrue(isinstance(util.getRootDataDirectory(), str))
 
  def testGetPath(self):
    if IGNORE_TEST:
      return
    SEP = os.path.join('x', 'y')[1]
    def test(path_initial, path_last):
      expected = ['a', 'b']
      if path_last is not None:
        expected.append(path_last)
      expected = SEP.join(expected)
      result = util.getPath(path_initial, path_last)
      self.assertTrue(all([x in result for x in path_initial]))
      self.assertEqual(result, expected)

    test(['a', 'b'], 'c')

  def testGetRootDataDirectory(self):
    if IGNORE_TEST:
      return
    result = util.getRootDataDirectory()
    split_result = result.split(DIR_SEPARATOR)
    self.assertTrue(
        len(util.DATA_DIRECTORIES.intersection(split_result)) > 0)

  def testGetDataPaths(self):
    if IGNORE_TEST:
      return
    path = util.getDataModelPath(cn.SQLDB_FILE)
    self.assertTrue(os.path.isfile(path))
    path = util.getDataModelPath(None)
    self.assertTrue(os.path.isdir(path))
    
  def testGetGeneNamesFromList(self):
    if IGNORE_TEST:
      return
    gene_column_list = ['MMP0335', 'DVU1485', 'DVU0168']
    column_list = copy.deepcopy(gene_column_list)
    column_list.append('rate')
    result_list = util.getGeneNamesFromList(column_list)
    self.assertEqual(set(gene_column_list),
        set(result_list))

  def testIsNanInDataFrame(self):
    if IGNORE_TEST:
      return
    length = 10
    df = pd.DataFrame({'a': list(range(length)), 
        'b': list(range(length))})
    b, columns = util.isNanInDataFrame(df)
    self.assertFalse(b)
    self.assertEqual(len(columns), 0)
    df['c'] = np.repeat(np.nan, length)
    b, columns = util.isNanInDataFrame(df)
    self.assertTrue(b)
    self.assertEqual(len(columns), 1)
    b, columns = util.isNanInDataFrame(df, nan_columns=['c'])
    self.assertFalse(b)
    self.assertEqual(len(columns), 0)

  def testGetColorMap(self):
    if IGNORE_TEST:
      return
    cmap = util.getColorMap()
    self.assertGreater(len(cmap), 0)

  def testGetColumnsFromDataFrame(self):
    if IGNORE_TEST:
      return
    data = list(range(10))
    col = 'a'
    df = pd.DataFrame({col: data, 'b': data})
    df2 = util.getColumnsFromDataFrame(df, [col])
    self.assertEqual(df2.columns.tolist(), [col])
    self.assertEqual(df2[col].tolist(), data)

  def testReplaceNan(self):
    if IGNORE_TEST:
      return
    col_a = 'a'
    col_b = 'b'
    length = 5
    def makeData(b_value):
      """
      Creates a dataframe and a deep copy with two columns
      :param object b_value:
      :return pd.DataFrame, pd.DataFrame:
      """
      a_data = range(length)
      b_data = np.repeat(b_value, length).tolist()
      df = pd.DataFrame({col_a: a_data, col_b: b_data})
      df_copy = copy.deepcopy(df)
      return df, df_copy

    df, df_copy = makeData(1) 
    util.replaceNan(df_copy)
    self.assertTrue(df.equals(df_copy))
    df, df_copy = makeData(np.nan,) 
    util.replaceNan(df_copy, value=2)
    self.assertEqual(df_copy[col_b][0], 2)
    df, df_copy = makeData(np.nan) 
    util.replaceNan(df_copy, columns=[col_a], value=2)
    self.assertTrue(np.isnan(df_copy[col_b][0]))

  def testChangeColumnValues(self):
    if IGNORE_TEST:
      return
    df = pd.DataFrame({'a': range(SIZE), 'b': range(SIZE)})
    dict_dc = {'a': 1, 'b': 2}
    func = lambda c, x: dict_dc[c]*x
    util.changeColumnValues(df, func)
    pairs = zip(df.a.tolist(), df.b.tolist())
    trues = [b == 2*a for a,b in pairs]
    self.assertTrue(all(trues))

  def testIsNull(self):
    if IGNORE_TEST:
      return
    self.assertTrue(util.isNull(np.nan))
    self.assertTrue(util.isNull(None))
    self.assertTrue(util.isNull("None"))
    self.assertTrue(util.isNull("none"))
    self.assertFalse(util.isNull("one"))
    self.assertFalse(util.isNull(range(10)))
    self.assertFalse(util.isNull(pd.DataFrame()))

  def testSetNoneList(self):
    if IGNORE_TEST:
      return
    self.assertEqual(len(util.setNoneList(None)), 0)
    self.assertEqual(len(util.setNoneList(range(10))), 10)

  def testGetColumnType(self):
    #if IGNORE_TEST:
    column_schemas = cn.TABLE_SCHEMAS.column_schemas.getSchemas()
    for column_schema in column_schemas:
      self.assertEqual(util.getColumnType(column_schema.name),
          column_schema.data_type)
    with self.assertRaises(RuntimeError):
      _ = util.getColumnType("DUMMY")

  def testStandardize(self):
    if IGNORE_TEST:
      return
    length = 10
    df = pd.DataFrame({
        'a': range(length),
        'b': np.repeat(1, length),
    })
    df_copy = copy.deepcopy(df)
    util.standardize(df, ['a', 'b'])
    self.assertTrue(np.isclose(df['a'].sum(), 0))
    self.assertTrue(np.isclose(df['b'].sum(), 0))
    self.assertTrue(np.isclose(np.std(df['a'].tolist()), 1))

  def testDeleteColumns(self):
    if IGNORE_TEST:
      return
    def test(columns):
      df = pd.DataFrame()
      for column in COLUMNS:
        df[column] = range(10)
      util.deleteColumns(df, columns)
      trues = [not x in df.columns for x in columns]
      self.assertTrue(all(trues))
      others = list(COLUMNS)
      [others.remove(x) for x in columns]
      trues = [x in df.columns for x in others]
      self.assertTrue(all(trues))

    test([COLUMNS[0]])
    test(COLUMNS[0:2])
    test(COLUMNS[1:3])

  def testTrimDF(self):
    if IGNORE_TEST:
      return
    columns = ['a', 'b', 'c', 'd']
    def test(keep_columns=None, delete_columns=None):
      df = pd.DataFrame()
      for column in columns:
        df[column] = range(10)
      df_new = util.trimDF(df, keep_columns=keep_columns,
          delete_columns=delete_columns)
      if keep_columns is not None:
        expecteds = keep_columns
      else:
        expecteds = list(columns)
        [expecteds.remove(c) for c in delete_columns]
      self.assertTrue(helpers.isValidDataFrame(df_new,
          expected_columns=expecteds))
      self.assertTrue(set(columns).issubset(df.columns))
      self.assertTrue(set(df.columns).issubset(columns))

    test(keep_columns=['a'])
    test(keep_columns=['a', 'b'])
    test(delete_columns=['b'])
    test(delete_columns=['b', 'a'])

  def testAddNullRow(self):
    if IGNORE_TEST:
      return
    df = pd.DataFrame({'a': range(SIZE), 'b': range(SIZE)})
    df1 = util.addNullRow(df)
    self.assertEqual(len(df1.index), SIZE+1)
    length = len([x for x in df1['a'] if util.isNull(x)])
    self.assertEqual(length, 1)

  def testGetDuplicates(self):
    if IGNORE_TEST:
      return
    self.assertEqual(util.getDuplicates([1,2, 1]), [1])
    self.assertEqual(util.getDuplicates([1,2]), [])
    self.assertEqual(util.getDuplicates([1,2, 1, 2]), [1, 2])

  def testTypeDF(self):
    if IGNORE_TEST:
      return
    Tester = collections.namedtuple('Tester',[
        'values', 'column', 'typ'])
    testers = [
        Tester(typ=bool, column=cn.IS_LOW_COVERAGE_ISOLATE,
             values=[1, True, None, np.nan]),
        Tester(typ=int, column=cn.POSITION,
             values=[1, 1.0, None, np.nan]),
        Tester(typ=str, column=cn.KEY_ISOLATE,
             values=['a', u'a', None, np.nan]),
        Tester(typ=float, column=cn.YIELD, 
             values=[64, 64.0, None, np.nan]),
        ]
    df = pd.DataFrame()
    for tester in testers:
      df[tester.column] = tester.values
    util.typeDF(df)
    for tester in testers:
      tester_values = tester.values
      df_values = df[tester.column].tolist()
      for nn in [0, 1]:
        if tester_values[nn] != df_values[nn]:
          import pdb; pdb.set_trace()
          pass
        self.assertEqual(tester_values[nn], 
            df_values[nn])
      for nn in [2, 3]:
        self.assertTrue(np.isnan(df_values[nn]))

  def testSelNonNull(self):
    if IGNORE_TEST:
      return
    values1 = [1, np.nan, 3, np.nan]
    values2 = [1, 2, np.nan, np.nan]
    self.assertEqual(util.selNonNull(values1, values2), 
        [1, 2, 3, np.nan])
    values1 = [10, np.nan, 3]
    values2 = [1, 2, np.nan]
    with self.assertRaises(ValueError):
      util.selNonNull([1, np.nan], [10, np.nan])

  def testMergeRowsColumns(self):
    if IGNORE_TEST:
      return
    common1 = [10, 20, 30, 40, 40]
    common2 = [10, 20, 30, 30, 50]
    only1 = [15, 25, 35, 45, 46]
    only2 = [115, 125, 135, 136, 145]
    merge_column1 = [1, 2, 3, 4, 4]
    merge_column2 = [1, 2, 3, 3, 5]
    df1 = pd.DataFrame({
        'merge_column': merge_column1,
        'common': common1,
        'only1': only1,
        })
    df2 = pd.DataFrame({
        'merge_column': merge_column2,
        'common': common2,
        'only2': only2,
        })
    df = util.mergeRowsColumns(df1, df2, 'merge_column')
    self.assertEqual(set(common1).union(common2),
        set(df['common']))
    self.assertEqual(set(merge_column1).union(merge_column2), 
        set(df['merge_column']))
    self.assertTrue(set(only1).issubset(df['only1']))
    self.assertTrue(set(only2).issubset(df['only2']))

  def testReadSQL(self):
    if IGNORE_TEST:
      return
    def test(cmd):
      df = util.readSQL(cmd)
      excludeds = [c for c in df.columns if c.count(':') > 0]
      self.assertEqual(len(excludeds), 0)
    #
    test("SELECT * FROM mutation;")
    test(cn.TABLE_SCHEMAS.getSchema(cn.TABLE_MUTATION))

  def testIsNumber(self):
    if IGNORE_TEST:
      return
    self.assertTrue(util.isNumber(1.0))
    self.assertTrue(util.isNumber('1.0'))
    self.assertTrue(util.isNumber(True))
    self.assertFalse(util.isNumber('a'))

  def testUnifyNullValuesList(self):
    if IGNORE_TEST:
      return
    def makeDF(values):
      return pd.DataFrame({'a': values})
    #
    values = range(10)
    df_value = makeDF(values)
    df_result = util.unifyNullValues(df_value)
    self.assertTrue(df_value.equals(df_result))
    df_value = pd.DataFrame(makeDF([0, np.nan, 2]))
    df_result = util.unifyNullValues(df_value)
    self.assertTrue(df_value.equals(df_result))
    values = [None, np.nan, 2.0]
    df_value = pd.DataFrame(makeDF(values))
    df_result = util.unifyNullValues(df_value)
    expected = list(values)
    expected[0] = cn.NONE
    df_expected = pd.DataFrame(makeDF(expected))
    self.assertTrue(df_result.equals(df_result))

  def testUnifyNullValuesDF(self):
    if IGNORE_TEST:
      return
    def testSeries(series1, series2):
      pairs = zip(series1.tolist(), series2.tolist())
      b = all([(x == y) or (np.isnan(x) and np.isnan(y))
               for x,y in pairs])
      self.assertTrue(b)

    df = pd.DataFrame({
        'a': range(3),
        'b': [0, np.nan, 2],
        'c': [None, np.nan, 2],
        })
    df_new = util.unifyNullValues(df)
    for column in ['a', 'b']:
      testSeries(df_new[column], df[column])
    expecteds = df['c'].tolist()
    expecteds[0] = cn.NONE
    df_new['expected'] = expecteds
    testSeries(df_new['c'], df_new['expected'])

  def testPruneNullRows(self):
    if IGNORE_TEST:
      return
    df = pd.DataFrame({
        'a': range(3),
        'b': range(3)
        })
    df_new = util.pruneNullRows(df)
    self.assertTrue(df.equals(df_new))
    df.loc[2,:] = cn.NONE
    df_new = util.pruneNullRows(df)
    df_expected = df.loc[0:1,:]
    self.assertTrue(df_expected.equals(df_new))

  def testPruneRowsWithNullColumns(self):
    if IGNORE_TEST:
      return
    df = pd.DataFrame({
        'a': range(3),
        'b': range(3)
        })
    df_new = util.pruneRowsWithNullColumns(df, ['a', 'b'])
    self.assertTrue(df.equals(df_new))
    df.loc[2, 'b'] = np.nan
    df_new = util.pruneRowsWithNullColumns(df, ['a', 'b'])
    self.assertEqual(len(df_new.index), len(df.index)-1)

  def testCoerceDF(self):
    if IGNORE_TEST:
      return
    self.df.loc[0, COL_A] = '0'
    df = util.coerceDF(self.df)
    self.assertEqual(df[COL_A].dtype, np.dtype('O'))
    # Rename to a known integer
    self.df.rename(columns={COL_A: cn.POSITION}, inplace=True)
    df = util.coerceDF(self.df)
    self.assertEqual(df[cn.POSITION].dtype, np.dtype('int64'))
    # Rename to a known bool
    self.df.rename(
        columns={cn.POSITION: cn.IS_AN_MUTATION}, inplace=True)
    self.df[cn.IS_AN_MUTATION] = True
    df = util.coerceDF(self.df)
    self.assertEqual(self.df[cn.IS_AN_MUTATION].dtype, np.dtype('bool'))

  def testAppendUnique(self):
    if IGNORE_TEST:
      return
    values = ['a', 'b']
    self.assertEqual(util.appendUnique(values, 'a'), values)
    values = ['a', 'b']
    new_values = util.appendUnique(values, 'c')
    self.assertTrue('c' in new_values)
    self.assertEqual(len(values) + 1, len(new_values))

  def testMergeLeft(self):
    if IGNORE_TEST:
      return
    def toString(a_list):
      return [str(x) for x in a_list]
  
    merge_column = 'b'
    df1 = pd.DataFrame({
      'a': [1, 2, 3],
      merge_column: [10, 20, np.nan],
      })
    df2 = pd.DataFrame({
      merge_column: [10, 30, np.nan],
      'c': [100, 200, 300],
      })
    df = util.mergeLeft(df1, df2, merge_column)
    self.assertEqual(toString(df1[merge_column]),
        toString(df[merge_column]))

  def testToNumber(self):
    if IGNORE_TEST:
      return
    self.assertEqual(util.toNumber('2'), 2)
    self.assertEqual(util.toNumber('2.0'), 2.0)
    self.assertTrue(np.isnan(util.toNumber('2.a')))

  def testUnpivot(self):
    if IGNORE_TEST:
      return
    DATA = range(5)
    LABEL_COLUMN = 'label_column'
    VALUE_COLUMN = 'value_column'
    df = pd.DataFrame({
        'a': DATA,
        'b': DATA,
        'c': DATA,
        })
    df_unpivot = util.unpivot(df, label_column=[LABEL_COLUMN],
        value_column=VALUE_COLUMN)
    self.assertTrue(helpers.isValidDataFrame(df_unpivot,
        [LABEL_COLUMN, VALUE_COLUMN]))

  def testMakeMatrix(self):
    if IGNORE_TEST:
      return
    df_base = util.readSQL("select * from genotype_phenotype")
    def test(row_name, column_name):
      df = df_base[[row_name, column_name]].copy()
      df[cn.COUNT] = 1
      df_matrix = util.makeMatrix(df, row_name=row_name,
          column_name=column_name)
      columns = list(df[column_name].unique())
      if np.nan in columns:
        columns.remove(np.nan)
      self.assertTrue(helpers.isValidDataFrame(df_matrix, columns))
    #
    test(cn.KEY_ISOLATE, cn.KEY_MUTATION)
    test(cn.LINE, cn.KEY_MUTATION)
    test(cn.KEY_ISOLATE, cn.GENE_ID)
    test(cn.LINE, cn.GENE_ID)
    test(cn.KEY_ISOLATE, cn.GENE_POSITION)
    test(cn.LINE, cn.GENE_POSITION)

  def testRemoveDuplicatesFromList(self):
    def test(values, trimmed_values):
      new_values = util.removeDuplicatesFromList(values)
      if isinstance(values[0], list):
        self.assertEqual(new_values, trimmed_values)
      else:
        self.assertEqual(set(new_values), set(trimmed_values))
    #
    test(['1', '2'], ['1', '2'])
    test([1, 2], [1, 2])
    test([1, 1, 2], [1, 2])
    test(['1', '1', '2'], ['1', '2'])
    test([['a', 'b'], ['c', 'd']],
        [['a', 'b'], ['c', 'd']])
    #
    test([['a', 'b'], ['c', 'd'], ['a', 'b']],
        [['a', 'b'], ['c', 'd']])

  def testNormalizedMean(self):
    VALUES = range(10)
    values = list(VALUES)
    values.append(100)
    mean3 = util.normalizedMean(values, 3)
    self.assertTrue(np.isclose(np.mean(VALUES), 
        np.std(VALUES)*mean3))
    mean1 = util.normalizedMean(values, 0.1)
    self.assertGreater(mean3, mean1)

  def testTrimExtremeValues(self):
    if IGNORE_TEST:
      return
    values = list(range(11))
    values.append(100)
    new_values1 = util.trimExtremeValues(values, 100)
    self.assertEqual(values, new_values1)
    new_values2 = util.trimExtremeValues(values, 3)
    self.assertTrue(set(new_values2).issubset(values))
    values = range(11)
    new_values3 = util.trimExtremeValues(values, 0)
    self.assertEqual(len(new_values3), 1)

  def testSelectRows(self):
    if IGNORE_TEST:
      return
    SIZE = 10
    DATA = range(SIZE)
    df = pd.DataFrame({
        COL_A: DATA,
        COL_B: DATA,
        })
    #
    df_new = df.copy(deep=True)
    df_new = util.selectRows(df_new, [])
    self.assertTrue(df_new.equals(df))
    #
    df_new = util.selectRows(df_new,
        [lambda r: r[COL_A] >= SIZE/2])
    self.assertEqual(len(df_new.index), SIZE/2) 

  def testFindRowColumn(self):
    if IGNORE_TEST:
      return
    size = 3
    df = pd.DataFrame({COL_A: range(size)})
    df.index = ['x', 'y', 'z']
    result = util.findRowColumn(df, lambda v: v < size)
    self.assertEqual(len(result), size)
    result = util.findRowColumn(df, lambda v: v < 0)
    self.assertEqual(len(result), 0)
    result = util.findRowColumn(df, lambda v: v < 1)
    self.assertEqual(len(result), 1)
    result = util.findRowColumn(df, lambda v: v < 2)
    self.assertEqual(len(result), 2)

  def testAggregateColumns(self):
    if IGNORE_TEST:
      return
    def aggregateFunc(df):
      # Initialization
      result = df[df.columns[0]]
      for col in df.columns:
        result = result + df[col]
      result = result.apply(lambda v: min(v, 1))
      return result
    #
    data = [1, 1, 0]
    df = pd.DataFrame({
      COL_A: data, 
      COL_B: data,
      COL_C: data,
      })
    new_column = util.aggregateColumns(df, [COL_A, COL_B], aggregateFunc)
    self.assertTrue(df[COL_C].equals(df[new_column]))

  def testFindCorrelatedColumns(self):
    if IGNORE_TEST:
      return
    data = range(5)
    df = pd.DataFrame({
        COL_A: data,
        COL_B: [2*d for d in data],
        COL_C: [1, -1, 1, -1, 0],
        })
    classes = util.findCorrelatedColumns(df)
    self.assertTrue(set([COL_A, COL_B]) in classes)
    self.assertEqual(len(classes), 1)

  def testSelectColumns(self):
    if IGNORE_TEST:
      return
    SET = set([COL_A, COL_B])
    constraint = lambda c: c in SET
    df = util.selectColumns(DATAFRAME, constraint)
    self.assertEqual(set(df.columns), SET)

  def testCorrelateWithPredictors(self):
    if IGNORE_TEST:
      return
    SIZE = 100
    makeData = lambda: [random.normalvariate(0, 1) 
        for _ in range(SIZE)]
    df_X = pd.DataFrame({
        COL_A: makeData(),
        COL_B: makeData(),
        COL_C: makeData(),
        })
    df_y = pd.DataFrame({
        cn.DEPVAR: [20*v + random.normalvariate(0,1) 
            for v in df_X[COL_A]],
        })
    df = util.correlateWithPredictors(df_X, df_y)
    self.assertEqual(len(df), len(df_X.columns))
    self.assertGreater(df.iloc[0, 0], 0.9)

  def testAppendWithColumnUnion(self):
    if IGNORE_TEST:
      return
    df1 = DATAFRAME.copy()
    df2 = DATAFRAME.copy()
    del df1[COL_C]
    del df2[COL_A]
    df = util.appendWithColumnUnion(df1, df2)
    self.assertTrue(helpers.isValidDataFrame(df,
        [COL_A, COL_B , COL_C]))
    self.assertEqual(len(df), 2*len(DATAFRAME))
    #
    df1 = pd.DataFrame()
    df = util.appendWithColumnUnion(df1, DATAFRAME)
    self.assertTrue(helpers.isValidDataFrame(df,
        [COL_A, COL_B , COL_C]))
    self.assertEqual(len(df), len(DATAFRAME))

  def testMakeProjectedPredictor(self):
    if IGNORE_TEST:
      return
    df_X = pd.DataFrame({
        COL_A: makeData(),
        COL_B: makeData(),
        })
    df_x = pd.DataFrame({
        COL_C: [2*x for x in df_X[COL_A]],
        })
    df = util.makeProjectedPredictor(df_X, df_x)
    self.assertTrue(helpers.isValidDataFrame(df,
        [COL_A, COL_B, COL_C]))
    for col in [COL_A, COL_B]:
        self.assertTrue(df_X[col].equals(df[col]))
    self.assertTrue(np.isclose(np.mean(df[COL_C]), 0.0))

  def testMakeProjectedPredictors(self):
    if IGNORE_TEST:
      return
    RTOL = 0.001
    isclose = lambda x,y: all([
        np.isclose(u, v, rtol=RTOL) for u, v in zip(x,y)])
    data = makeData()
    df_X = pd.DataFrame({
        COL_A: data,
        COL_B: makeData(),
        COL_C: [2*x for x in data],
        })
    df = util.makeProjectedPredictors(df_X)
    self.assertTrue(helpers.isValidDataFrame(df,
        [COL_A, COL_B, COL_C]))
    self.assertTrue(np.isclose(np.mean(df[COL_C]), 0.0, rtol=RTOL))

  def testExtendBin(self):
    if IGNORE_TEST:
      return
    result = util.extendBin(bin(2), 3)
    self.assertEqual(result, ['0', '1', '0'])
    result = util.extendBin(bin(2), 2)
    self.assertEqual(result, ['1', '0'])

  def testAddBinColumn(self):
    if IGNORE_TEST:
      return
    df = pd.DataFrame({
        COL_A: [0, 0, 1, 1],
        COL_B: [0, 1, 0, 1],
        })
    util.addBinColumn(df)
    self.assertTrue(helpers.isValidDataFrame(df,
        [COL_A, COL_B, cn.VALUE]))
    self.assertEqual(df[cn.VALUE].tolist(), df.index.tolist())

  def testnCr(self):
    if IGNORE_TEST:
      return
    self.assertEqual(util.nCr(5, 3), 10)
    self.assertEqual(util.nCr(5, 5), 1)

  def testFindKey(self):
    if IGNORE_TEST:
      return
    df = pd.DataFrame({
        COL_A: [0, 1, 0, 1, 0, 1, 0, 1],
        COL_B: [0, 0, 1, 1, 0, 0, 1, 1],
        COL_C: [0, 0, 0, 0, 1, 1, 1, 1],
        COL_D: [0, 1, 0, 1, 0, 1, 0, 1],
        })
    keys = util.findKey(df)
    self.assertEqual(keys, [COL_A, COL_B, COL_C])
    keys = util.findKey(df, required_columns=[COL_D])
    self.assertEqual(keys, [COL_B, COL_C, COL_D])

  def testMakeGeneDescriptionDF(self):
    if IGNORE_TEST:
      return
    ggene_ids = ['DVU0804', 'DVU1453', 'DVU2659', 'DVU2941',
        'DVH__.2982788.GENOME_CORRELATIONCC']
    df = util.makeGeneDescriptionDF(ggene_ids)
    self.assertTrue(helpers.isValidDataFrame(df,
      [cn.GGENE_ID, cn.GENE_DESC]))
    self.assertEqual(len(ggene_ids), len(df.index) + 1)

  def testGetRootDirectory(self):
    if IGNORE_TEST:
      return
    with self.assertRaises(ValueError):
      path = util.getIdentifiedDirectory(key_directory="dummy")
    path = util.getIdentifiedDirectory()
    self.assertGreater(len(path), 1)

  def testGetDBPath(self):
    if IGNORE_TEST:
      return
    path = util.getDBPath()
    self.assertTrue("microbepy" in path)
    

if __name__ == '__main__':
    unittest.main()
