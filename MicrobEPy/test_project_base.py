import project_base

import sys
import unittest

IGNORE_TEST = False
    

class TestFunctions(unittest.TestCase):

  def testGetProjectDirectory(self):
    path = project_base.getProjectDirectory()
    self.assertGreater(len(path), 1)

  def testAdPythonPaths(self):
    old_length = len(sys.path)
    project_base.addPythonPaths()
    self.assertGreater(len(sys.path), old_length)


if __name__ == '__main__':
    unittest.main()
