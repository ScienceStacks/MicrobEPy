"""Classification tree for binary valued variables."""


import __init__
import constants as cn
import binary_tree_model
import util

import copy
import numpy as np
import pandas as pd


# Hyperparameters
# Minimum accuracy by doing a split
MIN_INCR_SCORE = 0  # Minimum amount by which accuracy must improve
MAX_NODE_SCORE = 0.99  # Score beyond which there is no splitting
NO_SPLIT_MAX_NODE_SCORE = "max_node_score"
NO_SPLIT_MIN_INCR_SCORE = "min_incr_score"


class BinaryTreeClassification(binary_tree_model.BinaryTreeModel):
  """
  Does binary classification. y values are either 0 or 1.
  Predictors are binary as well.
  All nodes (except the root) classify either 0's or 1's
  for the depdendent variable, which is trained using df_y.
  This classification is recursively refined by adding
  columns from df_X, which is referred to as the split column. 
  self._side indicates whether the current
  node corresponds to the 0 or 1 values of the parent's
  split column. The state split value specifies that the
  node classifies values 0 in the dependent variable or 1.
  """

  def __init__(self, min_incr_score=MIN_INCR_SCORE, 
      max_node_score=MAX_NODE_SCORE, **kwargs):
    """
    :param float min_incr_score: minimal increase in a score
                                 to keep a subtree
    :param float max_node_score: Maximum node score for a split
    :param dict kwargs: optional arguments passed to parent
    """
    super(BinaryTreeClassification, self).__init__(**kwargs)
    self._min_incr_score = min_incr_score
    #
    self._max_node_score = max_node_score
    self.classify_value = np.nan

  def _printStats(self):
    return "ACC=%2.3f CLS=%f" % (self.node_score, 
        self.classify_value)

  @staticmethod
  def _calcSplitAccuracy(df_X, df_y):
    """
    Calculates the fraction of correct classifications (accuracy)
    for splitting on each column in df_X
    :param pd.DataFrame df_X:
    :param pd.DataFrame df_y:
    :return pd.Series: row index is column of df_X
    """
    # Calculate counts of correct classification
    width = len(df_X.columns)
    df_y_wide = pd.concat([df_y]*width, axis=1)
    df_y_wide.columns = df_X.columns
    df_both_one = df_X * df_y_wide
    #
    df_X_not = df_X.applymap(lambda v: 1 - v)
    df_y_wide_not = df_y_wide.applymap(lambda v: 1 - v)
    df_both_zero = df_X_not * df_y_wide_not
    ser_total_same = df_both_one.sum() + df_both_zero.sum()
    # Calculate classification accuracy
    length = len(df_X)
    ser_result = ser_total_same / length
    return ser_result

  def getScore(self):
    return self.score(None, None)

  def _findSplitColumn(self):
    """
    :return str, float, is_classify_one:
       str - split column, fraction correct
       float - split_score
    """
    cls = self.__class__
    # Find the column that maximizes the fraction correct classifications
    df_X_c = self.df_X.applymap(lambda v: 1 - v)
    # Compute the accuracy of predicting the dependent variable
    # Case 1: predict that y==1 if x==1; y==0 if x==0
    ser_q1 = cls._calcSplitAccuracy(self.df_X, self.df_y)
    # Case 0: predict that y==1 if x==0; y==0 if x==1
    ser_q0 = cls._calcSplitAccuracy(df_X_c, self.df_y)
    #
    ser_q = pd.Series([max(q0, q1) for q0, q1 in zip(ser_q0, ser_q1)])
    ser_q.index = ser_q0.index
    ser_q_sort = ser_q.sort_values()
    split_score = ser_q_sort.tolist()[-1]  # Fractional accuracy
    split_column = ser_q_sort.index[-1]
    #
    return split_column, split_score

  def fit(self, df_X, df_y, is_first_call=True):
    """
    Recursively find splits to form tree.
    Updates instance variables as appropriate.
    Does pruning of the tree after it is built.
    :param pd.DataFrame df_X: columns with binary values
    :param pd.DataFrame df_y: single column with float
    :param bool is_first_call: indicates first call to fit
    :return BinaryTreeClassification: fitted model
    Updates self.split_column, self.split_score, self.trees
    """
    cls = self.__class__
    self.validate()
    if isinstance(df_y, pd.Series):
      df_y = pd.DataFrame(pd.Series)
    col_y = df_y.columns[0]
    # Update state
    calcNodeScore = lambda: sum(self.df_y[self.df_y.columns[0]])   \
        / (1.0*len(self.df_y))
    self._updateState(df_X, df_y, calcNodeScore)
    length = len(self.df_y)
    one_count = self.df_y.sum().values[0]
    if one_count >= len(self.df_y)/2.0:
      self.classify_value = cn.ONE
      node_accuracy_numerator = one_count
    else:
      self.classify_value = cn.ZERO
      node_accuracy_numerator = length - one_count
    self.node_score = node_accuracy_numerator / (1.0*length)
    if self.node_score > self._max_node_score:
      # Don't split any further
      self._makeLeafState()
      self.no_split = NO_SPLIT_MAX_NODE_SCORE
    #
    columns = self.getSplitColumns()
    if len(columns) == 0:
      self.no_split = binary_tree_model.NO_SPLIT_NO_COLUMNS
      return copy.deepcopy(self)
    self.df_X = self.df_X[columns].copy()
    #
    split_column, split_score = self._findSplitColumn()
    # Acceptable split
    self.split_column = split_column
    self.split_score = split_score
    # Recursively do splits
    self._doSplits()
    # Prune the tree
    if is_first_call:
      self._prune()
    self.validate()
    return copy.deepcopy(self)

  def _prune(self):
    """
    Prunes the Tree based on whether the accuracy of
    an interior node is larger than the weighted accuracy
    of the leaves of the subtree for that interior node.
    """
    # Find the nodes to keep
    # The leaf should have the maximum score, and
    # this score should be sufficiently large.
    self.validate()
    all_keeps = []  # Keep the root
    leaves = self.findLeaves()
    for leaf in leaves:
      path = leaf.findPathToRoot()
      cur_leaf = leaf
      keeps = [cur_leaf]
      for node in path:
        incr_score =  node.getScore() - node.node_score
        if incr_score >= self._min_incr_score:
          keeps.append(node)
        else:
          cur_leaf = node
          keeps = [cur_leaf]
      all_keeps.extend(keeps)
    all_keeps = list(set(all_keeps))  # Eliminate duplicates
    # Delete nodes that are not to be kept
    remove_nodes = set(self.findNodes()).difference(all_keeps)
    for node in remove_nodes:
      if node in self.findNodes():
        parent = node.parent
        if parent is not None:
          parent.deleteSubtree()
          parent.no_split = NO_SPLIT_MIN_INCR_SCORE
    self.validate()
    return

  def score(self, _1, _2):
    """
    Computes the accuracy for this classification
    :param object _: placeholder signature is
                     compatible with scikit
    :return float:
    """
    score = 0.0
    leaves = self.findLeaves()
    for leaf in leaves:
      score += leaf.count*leaf.node_score
    score = score/self.count
    return score

  def getParameterDF(self):
    """
    Constructs a DataFrame with columnes for the variable name
    :return pd.DataFrame:
    """
    COL = 'dummy_column'
    df = pd.DataFrame({COL: [0]})
    nodes = set(self.findNodes())
    non_leaves = {n.split_column: n for n in nodes 
                 if n.split_column != cn.LEAF_NODE}
    for col in self.df_X:
      if col in non_leaves.keys():
        df[col] = 1
      else:
        df[col] = 0
    del df[COL]
    return df

  def calcPredictValue(self):
    """
    Provides the predicted value for the node.
    :return float:
    """
    return self.classify_value
