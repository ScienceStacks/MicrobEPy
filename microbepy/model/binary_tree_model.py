"""Abstract class for BinaryTreeRegression and BinaryTreeClassification"""

from microbepy.common import constants as cn
from microbepy.common import util

import copy
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score

SPACES = " "
BIN2STG = {cn.ZERO: "-", cn.ONE: "+", None: ""}
# Hyperparameters
MIN_SAMPLES = 3
MAX_DEPTH = 5  # Maximum depth of the tree
# Reasons for failing to split
NO_SPLIT_NO_COLUMNS = "no_columns"


class BinaryTreeModel(object):
  """
  A BinaryTreeModel is a tree where each node is itself a 
  BinaryTreeModel. A node corresponds to a column in the
  predictor matrix df_X, and each column is binary valued.
  """

  def __init__(self, min_samples=MIN_SAMPLES,
      max_depth=MAX_DEPTH, 
      cur_depth=0, 
      side=None,
      is_validate=False):
    """
    :param int min_samples: Minimum number of samples in a node
    :param int max_depth: How far the tree is expanded
    :param int cur_depth: Depth of this subtree in the tree
    :param int side: 0 or 1, reflecting whether the node classifies
                     0's or 1's
    """
    self._min_samples = min_samples
    self._max_depth = max_depth
    self._cur_depth = cur_depth
    self._side = side
    self.no_split = None  # Reason why node isn't split if not None
    self._is_validate = is_validate
    # State for each node
    self.parent = None
    self.df_X = None
    self.df_y = None
    self.count = 0  # Number of observations at this node
    self.node_score = np.nan  # Score before split
    self.mean = np.nan
    self.std = np.nan
    self._makeLeafState()

  def _makeLeafState(self):
    """
    Initializes state to values for a leaf node.
    """
    self.trees = {b: None for b in cn.BINARY_VALUES}
    self.split_score = np.nan  # Score after split
    self.incremental_score = np.nan  # Incremental score for split
    self.split_column = cn.LEAF_NODE  # Assume leaf node

  ##############################################################
  # Methods used by sub-classes: Graph manipulation
  ##############################################################
  def validate(self):
    """
    Validates the integrity of the tree.
    """
    if self._is_validate:
      for subtree in self.trees.values():
        if subtree is not None:
          if subtree.parent != self:
            print("Child does not point to parent!")
            import pdb; pdb.set_trace()
    bools = [t is None for t in self.trees.values()]
    if all(bools):
      if self.split_column != cn.LEAF_NODE:
        print("Leaf is improperly named.")
        import pdb; pdb.set_trace()
    elif not(any(bools)):
      if self.split_column == cn.LEAF_NODE:
        print("Non-leaf is improperly named.")
        import pdb; pdb.set_trace()

  def __str__(self):
    def makeSubtree(tree):
      if tree is None:
        return ""
      tree = "\n%s" % tree.__str__()
      return tree
    #
    if self._cur_depth == 0:
      result = "Total MSE: %f" % self.node_score
    else:
      result = ""
    spaces = ''.join([SPACES for _ in range(self._cur_depth)])
    stg = "%s%s %s CNT=%d" % (
        spaces, BIN2STG[self._side], self.split_column, self.count)
    full_stg = stg + " " + self._printStats()
    result = "%s\n%s" % (result, full_stg)
    for val in cn.BINARY_VALUES:
      result += "%s" % makeSubtree(self.trees[val])
    return result

  def findLeaves(self):
    """
    :return list-BinaryTreeRegression:
    """
    if self.split_column == cn.LEAF_NODE:
      result = [self]
    else:
      result = []
      for tree in self.trees.values():
        if tree is not None:
          result.extend(tree.findLeaves())
    return result

  def findNode(self, ser_x):
    """
    Finds the node that has the best match with ser_x values.
    :param pd.Series ser_x: indexed by split_columns, binary value
    :return BinaryTreeRegression:
    """
    col = self.split_column  # Split column for this node
    if col == cn.LEAF_NODE:
      return self
    val = ser_x[col]  # Which branch is used by ser_x
    if val not in cn.BINARY_VALUES:
      raise ValueError("Invalid binary value for %s: %f"
          % (col, val))
    next_tree = self.trees[val]
    if next_tree is None:
      return self
    else:
      return next_tree.findNode(ser_x)

  def findNodes(self, nodes=None, is_onlyleaves=False):
    """
    Finds all nodes in the Tree.
    :param list-BinaryTreeRegression nodes:
    :return list-BinaryTreeRegression:
    """
    nodes = util.setNoneList(nodes)
    if self.split_column == cn.LEAF_NODE:
      nodes.append(self)
      return nodes
    if not is_onlyleaves:
      nodes.append(self)
    for tree in self.trees.values():
      if tree is not None:
        tree.findNodes(nodes=nodes, is_onlyleaves=is_onlyleaves)
    return nodes
  
  def isRoot(self):
    return self.parent is None
  
  def isLeaf(self):
    return self.split_column == cn.LEAF_NODE

  def findPathToRoot(self):
    """
    :return list-BinaryTreeRegression:
    """
    cur_node = self
    nodes = [self]
    while not cur_node.isRoot():
      cur_node = cur_node.parent
      nodes.append(cur_node)
    return nodes

  def deleteChild(self, child):
    """
    Updates the state of the current node (self) when one
    of its children is deleted.
    :param BinaryTreeModel child:
    :return bool: True if found and removed
    """
    self.validate()
    found = False
    for k,v in self.trees.items():
      if v == child:
        found = True
        self.trees[k] = None
    if all([c is None for c in self.trees.values()]):
      self._makeLeafState()
    self.validate()
    return found

  def deleteSubtree(self, is_first_call=True):
    """
    Deletes all nodes that are children of this node, and makes this
    node into a root.
    :param bool is_first_call: Special handling for the root of subtree
    """
    self.validate()
    if is_first_call:
      # Transform the root of the subtree into a leaf
      self._makeLeafState()
    else:
      # Delete parent's connection with the child
      if self.parent is not None:
        found = parent.deleteChild(self)
        if not found:
          raise RuntimeError("Could not find child in parent trees.")
        self.parent = None  # Delete child connection with parent
    for tree in self.trees.values():
      if tree is not None:
        tree.deleteSubtree(is_first_call=False)
    self.validate()

  ##############################################################
  # Methods used by sub-classes: Model processing
  ##############################################################
  def getSplitColumns(self):
    """
    Finds the columns that are acceptable for splitting.
    :return list-str:
    """
    columns = []
    if self._cur_depth <= self._max_depth:
      # Check for too few samples remaining
      df_count = self.df_X.sum()
      max_samples = len(self.df_X) - self._min_samples
      columns = [c for c in df_count.index if
                (df_count[c] <= max_samples) and
                (df_count[c] >= self._min_samples)]
    return columns

  def _doSplits(self, **kwargs):
    """
    Split observations by doing fits recursively.
    """
    for val in cn.BINARY_VALUES:
      if self.split_column == cn.LEAF_NODE:
        break
      sel = [v == val for v in self.df_X[self.split_column]]
      is_single_value = all(sel) or (not any(sel))
      if not is_single_value:
        df_X_sub = self.df_X.loc[sel, :].copy()
        del df_X_sub[self.split_column]
        df_y_sub = self.df_y.loc[sel, :].copy()
        btr = self.__class__(min_samples=self._min_samples,
          max_depth=self._max_depth, cur_depth=self._cur_depth+1, 
          side=val, **kwargs)
        btr.fit(df_X_sub, df_y_sub, is_first_call=False)
        self.trees[val] = btr
        btr.parent = self
      else:
        self._makeLeafState()
    self.validate()

  def _updateState(self, df_X, df_y, calcNodeScore):
    """
    Updates state for the current node
    :param func calcNodeScore:
    """
    self.validate()
    cls = self.__class__
    if isinstance(df_y, pd.Series):
      df_y = pd.DataFrame(pd.Series)
    col_y = df_y.columns[0]
    #
    self.count = len(df_y)
    self.df_y = df_y.copy()
    self.df_X = df_X
    self.mean = self.df_y[col_y].mean()
    self.std = self.df_y[col_y].std()
    self.node_score = calcNodeScore()
    self.validate()

  def predict(self, df_X):
    """
    Constructs a prediction for a fitted model.
    :param pd.DataFrame df_X: binary predictors
    :return pd.DataFrame: cn.AVG, cn.STD, indexed as df_X
    """
    dfs = []
    for idx in range(len(df_X)):
      ser_x = df_X.iloc[idx, :]
      node = self.findNode(ser_x)
      df = pd.DataFrame({
          cn.AVG: [node.calcPredictValue()],
          cn.STD: [node.std],
          })
      dfs.append(df)
    df_result = pd.concat(dfs, sort=True)
    df_result.index = df_X.index
    return df_result

  def predictDF(self, df_X, df_y):
    """
    Constructs a dataframe to evaluate a model predictions.
    :param pd.DataFrame df_X:
    :param pd.DataFrame df_y:
    :return pd.DataFrame: cn.ESTIMATE, cn.OBSERVED, cn.RESIDUAL
    :raises ValueError:
    """
    ser_predicts = self.predict(df_X)[cn.AVG]
    col_y = df_y.columns[0]
    ser_observeds = df_y[col_y]
    ser_residuals = ser_observeds - ser_predicts
    df_predict = pd.DataFrame({
        cn.OBSERVED: ser_observeds,
        cn.ESTIMATE: ser_predicts,
        cn.RESIDUAL: ser_residuals,
        })
    return df_predict

  def getPredictors(self):
    """
    :return list-str: Predictor variables in the tree
    """
    nodes = self.findNodes()
    nodes = list(set(nodes))
    names = [n.split_column for n in nodes]
    names = [n for n in names if n != cn.LEAF_NODE]
    return names
