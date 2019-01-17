from microbepy.common import constants as cn
from microbepy.common import config
from microbepy.common import helpers
from microbepy.common import util

import os
import shutil
import unittest


IGNORE_TESTS = False


class TestFunctions(unittest.TestCase):

  def tearDown(self):
    if os.path.isdir(cn.CONFIG_DIR_PATH):
      shutil.rmtree(cn.CONFIG_DIR_PATH)

  def testInitialize(self):
    if IGNORE_TESTS:
      return
    config.initialize()
    self.assertTrue(os.path.isdir(cn.CONFIG_DIR_PATH))

  def testSetup(self):
    YAML_SETTING = {'x': 'dummy'}
    config.setup()
    self.assertTrue(os.path.isfile(cn.CONFIG_FILE_PATH))
    _ = config.setup(yaml_default=YAML_SETTING, is_forced=True)
    result = config.setup()
    self.assertEqual(result, YAML_SETTING)
    


if __name__ == '__main__': unittest.main()
