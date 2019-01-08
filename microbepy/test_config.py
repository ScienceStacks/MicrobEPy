import microbepy_init
import config

import os
import unittest

IGNORE_TEST = False
TEST_CONFIG_FILE = "config_file_test.yaml"
    

class TestFunctions(unittest.TestCase):

  def tearDown(self):
    if os.path.isfile(TEST_CONFIG_FILE):
      os.remove(TEST_CONFIG_FILE)

  def testGet(self):
    configuration = config.get()
    self.assertTrue("SQLDB_PATH" in configuration.keys())

  def testSet(self):
    KEY = "DUMMY"
    VALUE = 3
    config.set(KEY, VALUE, config_file=TEST_CONFIG_FILE)
    configuration = config.get(KEY, config_file=TEST_CONFIG_FILE)
    self.assertEqual(configuration, VALUE)
    
    


if __name__ == '__main__':
    unittest.main()
