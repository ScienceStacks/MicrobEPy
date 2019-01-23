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
CONFIG_DIR_PATH_SAV = "%s.sav" % cn.CONFIG_DIR_PATH


class TestFunctions(unittest.TestCase):

  def setUp(self):
    if os.path.isdir(CONFIG_DIR_PATH_SAV):
      shutil.rmtree(CONFIG_DIR_PATH_SAV)
    if os.path.isdir(cn.CONFIG_DIR_PATH):
      shutil.move(cn.CONFIG_DIR_PATH, CONFIG_DIR_PATH_SAV)

  def tearDown(self):
    if os.path.isdir(cn.CONFIG_DIR_PATH):
      shutil.rmtree(cn.CONFIG_DIR_PATH)
    if os.path.isdir(CONFIG_DIR_PATH_SAV):
      shutil.move(CONFIG_DIR_PATH_SAV, cn.CONFIG_DIR_PATH)

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
