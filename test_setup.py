import setup

import setuptools
import shutil
import sys
import unittest
from unittest import mock

IGNORE_TEST = False
    

class TestFunctions(unittest.TestCase):

  @mock.patch('subprocess.call', return_value=mock.Mock())
  def testDoRequirements(self, _):
    # Smoke test
    setup.pipInstall()

  @mock.patch('subprocess.call', return_value=mock.Mock())
  def testDoRequirements(self, _):
    # Smoke test
    setup.condaInstall()

  @mock.patch('os.mkdir', return_value=mock.Mock())
  @mock.patch('shutil.copyfile', return_value=mock.Mock())
  @mock.patch('subprocess.call', return_value=mock.Mock())
  @mock.patch('setuptools.setup', return_value=mock.Mock())
  def testMain(self, _, __, ___, ____):
    # Smoke test
    setup.main()


if __name__ == '__main__':
  unittest.main()
