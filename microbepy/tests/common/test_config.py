from microbepy.common import constants as cn
from microbepy.common import config
from microbepy.common import helpers
from microbepy.common import util

import builtins
import os
import shutil
import yaml
import unittest
from unittest import mock


IGNORE_TESTS = False


class TestFunctions(unittest.TestCase):

  def tearDown(self):
    if os.path.isdir(cn.CONFIG_DIR_PATH):
      shutil.rmtree(cn.CONFIG_DIR_PATH)

  @mock.patch('os.mkdir', return_value=mock.Mock())
  def testInitialize(self, _):
    if IGNORE_TESTS:
      return
    config.initialize()

  @mock.patch('os.mkdir', return_value=mock.Mock())
  @mock.patch('os.path.isfile', return_value=False)
  @mock.patch('yaml.dump', return_value=mock.Mock())
  def testSetup(self, _, __, ___):
    if IGNORE_TESTS:
      return
    YAML_SETTING = {'x': 'dummy'}
    config.setup()
    # self.assertTrue(os.path.isfile(cn.CONFIG_FILE_PATH))
    # _ = config.setup(yaml_default=YAML_SETTING, is_forced=True)
    # result = config.setup()
    # self.assertEqual(result, YAML_SETTING)
    


if __name__ == '__main__': unittest.main()
