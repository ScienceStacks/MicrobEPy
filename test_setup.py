import setup
import distutils.core

import shutil
import sys
import unittest
from unittest import mock

IGNORE_TEST = False
    

class TestFunctions(unittest.TestCase):

  @mock.patch('os.mkdir', return_value=mock.Mock())
  @mock.patch('shutil.copyfile', return_value=mock.Mock())
  def testCopyData(self, _, __):
    # Smoke test
    setup.copyData()

  @mock.patch('subprocess.call', return_value=mock.Mock())
  def testDoRequirements(self, _):
    # Smoke test
    setup.doRequirements()

  @mock.patch('os.mkdir', return_value=mock.Mock())
  @mock.patch('shutil.copyfile', return_value=mock.Mock())
  @mock.patch('subprocess.call', return_value=mock.Mock())
  @mock.patch('distutils.core.setup', return_value=mock.Mock())
  def testMain(self, _, __, ___, ____):
    # Smoke test
    setup.main()


if __name__ == '__main__':
  unittest.main()
