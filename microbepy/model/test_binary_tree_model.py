import microbepy_init
import constants as cn
import helpers
import util
import binary_tree_regression as btr
from helpers_test_model import makeTree
import binary_tree_model as btm
import util_data as ud

import copy
import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
TREE_SIZE = 5
COL_A = 'a'
COL_B = 'b'
COL_C = 'c'
SIZE_A = 3
SIZE_B = 2
DF_X = pd.DataFrame({
    COL_A: [1, 1, 1, 0],
    COL_B: [1, 1, 0, 0],
    })
DF_Y = pd.DataFrame({
    cn.VALUE: range(4)
    })
SMALL = 0.33
LARGE = 0.66
NAMES = ['a', 'b', 'c', 'd', 'e', 'f']
DEPTH = len(NAMES) - 1



################## TEST FUNCTIONS #################
class TestBinaryTreeModel(unittest.TestCase):

  def setUp(self):
    self.bt_regress = btr.BinaryTreeRegression(min_samples=0,
        is_validate=True)

  def testStr(self):
    if IGNORE_TEST:
      return
    tree = makeTree()
    tree_str = tree.__str__()
    self.assertEqual(tree_str.count(" a "), 1)
    for depth in range(1, DEPTH-1):
      count = 2**(depth-1)
      stg = "%s1 " % NAMES[depth]
      self.assertEqual(tree_str.count(stg), count)

  def testFindNode(self):
    if IGNORE_TEST:
      return
    def test(val):
      index = [NAMES[0]]
      for name in NAMES[1:]:
        new = index[-1] + name + str(val)
        index.append(new)
      ser_x = pd.Series([val for _ in range(len(NAMES))])
      ser_x.index = index
      #
      tree = makeTree()
      node = tree.findNode(ser_x)
      self.assertEqual(node.split_column, cn.LEAF_NODE)
    #
    test(0)
    test(1)

  def testPredict(self):
    if IGNORE_TEST:
      return
    self.bt_regress.fit(DF_X, DF_Y)
    df = self.bt_regress.predict(DF_X)
    self.assertTrue(helpers.isValidDataFrame(df,
        [cn.AVG, cn.STD], nan_columns=[cn.STD]))
    self.assertEqual(len(df), len(DF_X))

  def testFindNodes(self):
    if IGNORE_TEST:
      return
    self.bt_regress.fit(DF_X, DF_Y)
    leaves = self.bt_regress.findNodes(is_onlyleaves=True)
    trues = [t.split_column == cn.LEAF_NODE for t in leaves]
    self.assertTrue(all(trues))
    self.assertEqual(len(leaves), 3)
    #
    nodes = self.bt_regress.findNodes(is_onlyleaves=False)
    columns = [t.split_column for t in nodes 
               if t.split_column != cn.LEAF_NODE]
    self.assertEqual(set(columns), set(DF_X.columns))
    self.assertEqual(len(nodes), 5)

  def testPredictDF(self):
    if IGNORE_TEST:
      return
    self.bt_regress.fit(DF_X, DF_Y)
    self.bt_regress.predictDF(DF_X, DF_Y)
    df = self.bt_regress.predictDF(DF_X, DF_Y)
    self.assertTrue(helpers.isValidDataFrame(df,
        [cn.OBSERVED, cn.ESTIMATE, cn.RESIDUAL]))

  def testFindLeaves(self):
    self.bt_regress.fit(DF_X, DF_Y)
    leaves = self.bt_regress.findLeaves()
    self.assertEqual(len(leaves), 3)
    self.bt_regress.validate()

  def testFindPathToRoot(self):
    if IGNORE_TEST:
      return
    tree = makeTree(max_depth=TREE_SIZE)
    path = tree.findPathToRoot()
    self.assertEqual(len(path), 1)
    #
    leaves = tree.findLeaves()
    leaf = leaves[0]
    path = leaf.findPathToRoot()
    self.assertEqual(len(path), TREE_SIZE+1)
    #
    self.assertEqual(leaf, path[0])
    self.assertEqual(tree, path[-1])
    #
    for idx, tree in enumerate(path[:-1]):
      self.assertEqual(tree.parent, path[idx+1]) 

  def testDeleteChild(self):
    if IGNORE_TEST:
      return
    def testLeaf(parent, val):
      leaf = parent.trees[val]
      found = parent.deleteChild(leaf)
      self.assertTrue(found)
      self.assertFalse(leaf in parent.trees.values())
    #
    tree = makeTree(max_depth=TREE_SIZE)
    leaf = tree.findLeaves()[0]
    parent = leaf.parent
    #
    testLeaf(parent, 0)
    self.assertNotEqual(parent.split_column, cn.LEAF_NODE)
    testLeaf(parent, 1)
    self.assertEqual(parent.split_column, cn.LEAF_NODE)

  def testDeleteSubtree(self):
    if IGNORE_TEST:
      return
    OFFSET = 2
    tree = makeTree(max_depth=TREE_SIZE)
    old_tree = copy.deepcopy(tree)
    all_names = set([n.split_column for n in tree.findNodes()])
    leaf = tree.findLeaves()[0]
    path = leaf.findPathToRoot()
    subtree = path[TREE_SIZE - OFFSET - 1]
    remove_names = ([n.split_column for n in subtree.findNodes()])
    subtree.deleteSubtree()
    for name in remove_names:
      if name != cn.LEAF_NODE:
        self.assertFalse(name in str(tree))
    expecteds = all_names.difference(remove_names)
    for name in expecteds:
      self.assertTrue(name in str(tree))

  def testGetPredictors(self):
    if IGNORE_TEST:
      return
    tree = makeTree(max_depth=TREE_SIZE)
    predictors = tree.getPredictors()
    tot = 0  # Number of interior nodes
    for n in range(TREE_SIZE):
      tot += 2**n
    self.assertEqual(len(predictors), tot)
    

if __name__ == '__main__':
    unittest.main()
