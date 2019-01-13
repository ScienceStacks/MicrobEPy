from microbepy.common import helpers
from microbepy.common import util
from microbepy.common import constants as cn
from microbepy.common.group_collection  \
    import Group, GroupCollection, ELEMENT_SEPARATOR

import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
DATA1 = ['a', 'b', 'c']
DATA2 = ['d', 'b']
DATA_LIST = [DATA1, DATA2]
ALL_DATA = set(DATA1).union(DATA2)
TWO = 2
GROUP1 = Group(DATA1)
GROUP2 = Group(DATA2)
PREFIX = "prefix"


class TestGroup(unittest.TestCase):

  def setUp(self):
    self.group = Group(DATA1)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(self.group.len(), len(DATA1))
  
  def testCopy(self):
    if IGNORE_TEST:
      return
    group = self.group.copy()
    self.assertTrue(self.group.equals(group))

  def testRepr(self):
    if IGNORE_TEST:
      return
    rep = str(self.group)
    self.assertTrue(ELEMENT_SEPARATOR in rep)

  def testIntersection(self):
    if IGNORE_TEST:
      return
    group = Group(DATA2)
    new_group = self.group.intersection(group)
    self.assertEqual(new_group.len(), 1)

  def testEquals(self):
    if IGNORE_TEST:
      return
    self.assertTrue(self.group.equals(self.group))
    new_group = Group(self.group.objects)
    self.group.value = 1.0
    self.group.label = "Dummy"
    self.assertTrue(self.group.equals(self.group))

  def testMakeGroupFromString(self):
    if IGNORE_TEST:
      return
    rep = str(self.group)
    group = Group.makeGroupFromString(rep)
    self.assertTrue(self.group.equals(self.group))

  def testUnion(self):
    if IGNORE_TEST:
      return
    group = self.group.union(Group(DATA2))
    self.assertEqual(group.len(), len(ALL_DATA))


class TestGroupCollection(unittest.TestCase):

  def setUp(self):
    self.group_collection = GroupCollection()

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(isinstance(self.group_collection,
        GroupCollection))
    self.assertEqual(self.group_collection.len(), 0)

  def testAdd(self):
    if IGNORE_TEST:
      return
    self.group_collection.add(DATA1)
    self.assertEqual(self.group_collection.len(), 1)

  def testIntersection(self):
    if IGNORE_TEST:
      return
    self.group_collection.add(DATA1)
    other = GroupCollection([Group(DATA2)])
    groups = self.group_collection.intersection(other)
    self.assertEqual(groups.len(), 1)

  def testFlatten(self):
    if IGNORE_TEST:
      return
    self.group_collection.add(DATA1)
    self.group_collection.add(DATA2)
    group = self.group_collection.flatten()
    self.assertEqual(group.len(), len(ALL_DATA))

  def setupGroupCollection(self, group_collection=None, 
      is_label=False):
    if group_collection is None:
      group_collection = self.group_collection
    group1 = Group(DATA1)
    group2 = Group(DATA2)
    group1.value = 1
    group2.value = 2
    self.group_collection.add(group1)
    self.group_collection.add(group2)
    if is_label:
      self.group_collection.setGroupLabels(prefix=PREFIX)

  def testMakeValueDF(self):
    if IGNORE_TEST:
      return
    self.setupGroupCollection()
    df = self.group_collection.makeValueDF()
    self.assertEqual(len(df.columns),
        len(self.group_collection.groups))
    self.assertEqual(len(df.index),
        self.group_collection.flatten().len())
    for group in self.group_collection.groups:
      values = df[str(group)]
      count = len([v for v in values if not np.isnan(v)])
      self.assertEqual(count, group.len())
    # Test force_value. Ensure get the same results as before.
    df_new = self.group_collection.makeValueDF(force_value=1.0)
    df_bool = df.applymap(lambda v: not np.isnan(v))
    df_new_bool = df_new.applymap(lambda v: np.isclose(v, 1.0))
    for column in df_bool.columns:
      trues = df_bool[column] == df_new_bool
      self.assertTrue(all(trues))

  def testUnionDisjoint(self):
    if IGNORE_TEST:
      return
    # Setup
    self.setupGroupCollection()
    #
    group1 = Group(DATA1)
    group2 = Group(DATA2)
    group1.value = 1
    group2.value = 2
    group_collection1 = GroupCollection()
    group_collection1.add(group1)
    group_collection2 = GroupCollection()
    group_collection2.add(group2)
    # Test
    group_collection = group_collection2.unionDisjoint(
        group_collection1)
    self.assertTrue(self.group_collection.equals(group_collection))
    self.assertFalse(self.group_collection.equals(group_collection1))
    with self.assertRaises(ValueError):
      group_collection = group_collection2.unionDisjoint(
          group_collection2)

  def testEquals(self):
    if IGNORE_TEST:
      return
    self.setupGroupCollection()
    self.assertTrue(self.group_collection.equals(
        self.group_collection))

  def testCopy(self):
    if IGNORE_TEST:
      return
    self.setupGroupCollection()
    groups = self.group_collection.copy()
    self.assertTrue(self.group_collection.equals(
        groups))

  def testSetGroupLabels(self):
    if IGNORE_TEST:
      return
    self.setupGroupCollection(is_label=True)
    trues = [PREFIX in g.label for g in self.group_collection.groups]
    self.assertTrue(trues)

  def testGetGroupLabels(self):
    self.setupGroupCollection(is_label=True)
    labels = self.group_collection.getGroupLabels()
    self.assertEqual(len(labels), len(self.group_collection.groups))
    

if __name__ == '__main__':
    unittest.main()
