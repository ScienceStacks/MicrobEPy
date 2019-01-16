"""Tests for Helpers"""

from microbepy.common import helpers as he

import numpy as np
import pandas as pd
import os
import unittest

IGNORE_TEST = False
he.DEBUG = False  # Don't use pdb when invalid dataframe


class TestFunctions(unittest.TestCase):

  def testisValidDataFrame(self):
    length = 10
    #
    data = {'a': range(length), 'b': range(length)}
    columns = list(data.keys())
    df = pd.DataFrame(data)
    self.assertTrue(he.isValidDataFrame(df, list(data.keys()),
        key=['a', 'b']))
    #
    data = {'a': range(length), 'b': range(length)}
    df = pd.DataFrame(data)
    self.assertTrue(he.isValidDataFrame(df, 
        expected_columns=columns))
    #
    df2 = pd.DataFrame(df['a'])
    self.assertFalse(he.isValidDataFrame(df2, 
        expected_columns=columns))
    #
    data.update({'c': np.repeat(np.nan, length).tolist()})
    columns = list(data.keys())
    df3 = pd.DataFrame(data)
    self.assertFalse(he.isValidDataFrame(df3, columns))
    

if __name__ == '__main__':
    unittest.main()

