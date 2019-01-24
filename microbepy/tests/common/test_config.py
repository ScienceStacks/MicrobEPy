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
KEY = "key"
VALUE = "dummy"
YAML_DICT = {KEY: VALUE}


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

  def testSetupAndGet(self):
    if IGNORE_TESTS:
      return
    config.setup(yaml_dict=YAML_DICT)
    self.assertTrue(os.path.isfile(cn.CONFIG_FILE_PATH))
    value =  config.get(KEY)
    self.assertTrue(VALUE, value)
    

  def testGet(self):
    if IGNORE_TESTS:
      return
    with self.assertRaises(KeyError):
      value =  config.get(KEY)

    


if __name__ == '__main__': unittest.main()
