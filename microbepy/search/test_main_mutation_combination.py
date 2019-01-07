import microbepy_init
import helpers
import util
import constants as cn
from main_mutation_combination import main, doCombination
import mutation_context

import numpy as np
import os
import pandas as pd
import shutil
import unittest


IGNORE_TEST = False
DEPVAR = cn.RATE
MUTATION_COLUMN = cn.POSITION
TEST_FILE_PREFIX = "testcsv"
TEST_FILE = TEST_FILE_PREFIX + "do_combination.csv"
DONE_PROB = 0.2
MUTATION_CONTEXT = mutation_context.MutationContext(
    DEPVAR, MUTATION_COLUMN)


############################################################
def deleteFiles():
  files = os.listdir(os.curdir)
  for ffile in files:
    if TEST_FILE_PREFIX in ffile:
      shutil.move(ffile, os.path.join("/tmp", ffile))


############################################################
# Mock for MutationCombination
############################################################
class MutationCombinationMock(object):

  def __init__(self, context):
    if isinstance(context, str):
      import pdb; pdb.set_trace()
    self._context = context

  @classmethod
  def isMock(cls):
    """
    Used for dependency injection.
    """
    return True

  @classmethod
  def getFilePrefix(cls):
    """
    Used for dependency injection.
    """
    return "testcsv"
  
  def do(self, max_combination, is_tstat=False, 
      is_resample=False, lines=None, excludes=None,
      **kwargs):
    combinations = [['a', 'b']]
    if excludes == combinations:
      df = pd.DataFrame()
    else:
      df = pd.DataFrame()
      for line in [cn.LINE_ALL, cn.LINE_HA2]:
        df = df.append(pd.DataFrame({
            cn.DEPVAR: [self._context.depvar],
            cn.MUTATION_COLUMN: [self._context.mutation_column],
            cn.MUTATIONS: [str(c) for c in combinations],
            cn.SL_TSTAT: [1.0],
            cn.SL_RESAMPLE: [1.0],
            cn.VALUE: [1.0],
            cn.COUNT: [1],
            cn.GROUP: ["0b10"],
            cn.LINE: [line],
            'KWARGS': kwargs.values(),
            }))
    return df, df
 


############################################################
# Tests for MutationCombination
############################################################
class TestMainMutationCombination(unittest.TestCase):

  def setUp(self):
    deleteFiles()

  def tearDown(self):
    deleteFiles()

  def testDoCombination(self):
    if IGNORE_TEST:
      return
    files = [TEST_FILE, TEST_FILE]
    def test(expected_lines):
      df_min, df_max =  \
          doCombination(1, MUTATION_CONTEXT, [cn.LINE_ALL],
          files, combination_class=MutationCombinationMock)
      for df in [df_min, df_max]:
        self.assertEqual(len(df), expected_lines)
      return df_min, df_max
    #
    df_min, df_max = test(2)
    #
    df_min.to_csv(TEST_FILE)
    df_min, df_max = test(0)

  def testMain(self):
    if IGNORE_TEST:
      return
    main(combination_class=MutationCombinationMock)
    files = os.listdir(os.curdir)
    files = [f for f in files if TEST_FILE_PREFIX in f]
    expected = 2*len(cn.DEPVARS)*len(cn.MUTATION_COLUMNS)
    self.assertEqual(len(files), expected)


if __name__ == '__main__':
    unittest.main()
