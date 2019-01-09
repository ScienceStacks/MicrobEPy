import setup

import shutil
import sys
import unittest
from unittest import mock

IGNORE_TEST = False
    

class TestFunctions(unittest.TestCase):

  @mock.patch('os.mkdir', return_value=mock.Mock())
  @mock.patch('shutil.copyfile', return_value=mock.Mock())
  def testCopyData(self, _, __):
    setup.copyData()


if __name__ == '__main__':
  unittest.main()
