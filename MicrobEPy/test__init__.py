import __init__

import sys
import unittest

IGNORE_TEST = False
    

class TestFunctions(unittest.TestCase):

  def testGetProjectDirectory(self):
    path = __init__.getProjectDirectory()
    self.assertGreater(len(path), 1)

  def testAdPythonPaths(self):
    old_length = len(sys.path)
    __init__.addPythonPaths()
    self.assertGreater(len(sys.path), old_length)

  def testPath(self):
    for path in sys.path:
      if __init__.PROJECT_NAME in path:
        self.assertTrue(path.count(__init__.PROJECT_NAME) in [1, 2])


if __name__ == '__main__':
    unittest.main()
