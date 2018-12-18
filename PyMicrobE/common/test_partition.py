import __init__
from partition import Partition

import copy
import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
ELEMENTS = ['a', 'b', 'c', 'd', 'e']
SET1 = set(['a', 'b', 'c'])
SET2 = set(['d', 'e'])
PARTITION = [SET1, SET2]


class TestPartition(unittest.TestCase):

  def setUp(self):
    self.partition = Partition(ELEMENTS)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(set(ELEMENTS), self.partition.elements)
    with self.assertRaises(ValueError):
      Partition(ELEMENTS, [])

  def testIsPresent(self):
    if IGNORE_TEST:
      return
    a_set = set(ELEMENTS[0])
    self.assertTrue(self.partition.isPresent(a_set))
    a_set = set("BAD")
    self.assertFalse(self.partition.isPresent(a_set))
    self.partition = Partition(ELEMENTS, partition=PARTITION)
    self.assertTrue(self.partition.isPresent(SET1))

  def testFindSet(self):
    self.partition = Partition(ELEMENTS, partition=PARTITION)
    idx1 = self.partition._findSet(SET1)
    self.assertTrue(idx1 in range(len(PARTITION)))
    idx2 = self.partition._findSet(SET2)
    self.assertTrue(idx2 in range(len(PARTITION)))
    self.assertTrue(idx1 != idx2)

  def testMove(self):
    element = ELEMENTS[0]
    self.partition = Partition(ELEMENTS, partition=PARTITION)
    set_mutation = self.partition.move(element, SET2,
        is_update=False)
    self.assertEqual(set_mutation.new_src, 
        SET1.difference(element))
    self.assertEqual(set_mutation.cur_dst, SET2)
    self.assertEqual(set_mutation.new_dst, 
        SET2.union(element))
    self.assertEqual(self.partition.sets, PARTITION)
    #
    set_mutation = self.partition.move(element, SET2,
        is_update=True)
    self.assertEqual(set_mutation.new_src, 
        SET1.difference(element))
    self.assertEqual(set_mutation.cur_dst, SET2)
    self.assertEqual(set_mutation.new_dst, 
        SET2.union(element))
    self.assertEqual(self.partition.sets,
        [set_mutation.new_src, set_mutation.new_dst])

  def testFindSetWithElement(self):
    self.partition = Partition(ELEMENTS, partition=PARTITION)
    a_set = self.partition.findSetWithElement(list(SET1)[0])
    self.assertEqual(a_set, SET1)
    with self.assertRaises(ValueError):
      a_set = self.partition.findSetWithElement('DUMMY',
          is_present=True)
    

if __name__ == '__main__':
    unittest.main()
