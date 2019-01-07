import microbepy_init

import sys
import unittest

IGNORE_TEST = False
    

class TestFunctions(unittest.TestCase):

  def testGetProjectDirectory(self):
    path = microbepy_init.getProjectDirectory()
    self.assertGreater(len(path), 1)

  def testAdPythonPaths(self):
    old_length = len(sys.path)
    microbepy_init.addPythonPaths()
    self.assertGreater(len(sys.path), old_length)

  def testPath(self):
    for path in sys.path:
      if microbepy_init.PROJECT_NAME in path:
        self.assertTrue(path.count(microbepy_init.PROJECT_NAME) in [1, 2])


if __name__ == '__main__':
    unittest.main()
