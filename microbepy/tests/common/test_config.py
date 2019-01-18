from microbepy.common import constants as cn
from microbepy.common import config
from microbepy.common import helpers
from microbepy.common import util

import builtins
import os
import shutil
import yaml
import unittest


IGNORE_TESTS = False


class TestFunctions(unittest.TestCase):

  def tearDown(self):
    if os.path.isdir(cn.CONFIG_DIR_PATH):
      config.setup(yaml_default=cn.YAML_DEFAULT, is_forced=True)

  def testInitialize(self):
    if IGNORE_TESTS:
      return
    config.initialize()

  def testSetup(self):
    if IGNORE_TESTS:
      return
    YAML_SETTING = {'x': 'dummy'}
    config.setup()
    if True:
      self.assertTrue(os.path.isfile(cn.CONFIG_FILE_PATH))
      _ = config.setup(yaml_default=YAML_SETTING, is_forced=True)
      result = config.setup()
      self.assertEqual(result, YAML_SETTING)
    


if __name__ == '__main__': unittest.main()
