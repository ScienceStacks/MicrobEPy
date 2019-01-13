from microbepy.common.equivalence_class import EquivalenceClass
import numpy as np
import unittest


IGNORE_TEST = False
SIZE = 20
ITERABLES = range(SIZE)
MODULO = 4
RELATION = lambda x,y: (x % MODULO) == (y % MODULO)
    

class TestGenomeRegression(unittest.TestCase):

  def setUp(self):
    self.equiv = EquivalenceClass(ITERABLES, RELATION)

  def testConstructor(self):
    self.assertEqual(len(self.equiv.classes), 0)

  def testDo(self):
    self.equiv.do()
    self.assertEqual(len(self.equiv.classes), MODULO)
    self.assertEqual(len(self.equiv.classes[0]), SIZE/MODULO)

  def testValidate(self):
    self.equiv.do()
    self.equiv.validate()
    self.equiv.classes[0] = [1, 2]
    with self.assertRaises(ValueError):
      self.equiv.validate()
    


if __name__ == '__main__':
    unittest.main()
