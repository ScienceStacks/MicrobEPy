"""Helpers for testing the models directory."""

import binary_tree_regression as btr
import binary_tree_model as btm
import constants as cn

import pandas as pd
import numpy as np

DEPTH = 4
NAMES = ['a', 'b', 'c', 'd', 'e', 'f']

################# HELPERS ######################
def makeTree(cls=btr.BinaryTreeRegression, max_depth=DEPTH):
  def setValues(tree, side):
    tree.df_y = pd.DataFrame(np.repeat(1, 10))
    tree._side = side
    tree.split_score = 1.0
    tree.node_score = 1.0
  def addTree(tree, side, depth):
    new_tree = cls()
    setValues(new_tree, side)
    new_tree._cur_depth = depth + 1
    if new_tree._cur_depth < max_depth:
      new_tree.split_column = tree.split_column + NAMES[depth+1]  \
          + str(side)
    return new_tree
  def buildTree(tree):
    if tree._cur_depth < max_depth:
      for val in cn.BINARY_VALUES:
        side = val
        tree.trees[val] = addTree(tree, side, tree._cur_depth)
        tree.trees[val].parent = tree
        buildTree(tree.trees[val])
  def assignCounts(tree):
    """
    Count is the number of leaves in the subtree
    """
    if tree.isLeaf():
      tree.count = 1
    tree.count = len(tree.findLeaves())
    for val in cn.BINARY_VALUES:
      sub_tree = tree.trees[val]
      if sub_tree is not None:
        assignCounts(sub_tree)
  #
  root = cls()
  setValues(root, cn.ONE)
  root.split_column = 'a'
  buildTree(root)
  assignCounts(root)
  return root
