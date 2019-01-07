"""
Describes Groups of Groups.

A Group is a set of objects that have unique string representations.
A Group has an unique ID or representation that is formed from the string representation of
its objects. A Group may also have a value and a label; these attributes
are used when plotting groups.
A GroupCollection is a set of Groups.
"""


import microbepy_init
import constants as cn
import util

from collections import namedtuple
import numpy as np
import pandas as pd

ELEMENT_SEPARATOR = "--"


class Group(object):
  """A group are similar objects nested in a collection."""

  def __init__(self, objects, value=np.nan, label=None):
    """
    :param collection-objects objects: must be iterables
    :param float value:
    :param str label:
    Note: objects must be unique
    """
    self.objects = [o for o in objects]
    self.objects.sort()
    if len(self.objects) > len(set(self.objects)):
      raise ValueError("Objects are not unique")
    self.value = value
    self.label = label

  def __repr__(self):
    label = [str(o) for o in self.objects]
    return ELEMENT_SEPARATOR.join(label)

  def id(self):
    return self.__repr__()

  def intersection(self, other):
    """
    Construct a group that is the intersection of another group
    Notes:
      1. Does not perserve value or label
    """
    return Group(set(self.objects).intersection(other.objects))

  @staticmethod
  def makeGroupFromString(group_string):
    """
    Constructs a group from a string that delinates members by the ELEMENT_SEPARATOR.
    :param str group_string:
    :return Group:
    """
    return Group(group_string.split(ELEMENT_SEPARATOR))

  def equals(self, other):
    """
    :param Group other:
    :return bool:
    """
    objects = set(self.objects)
    if len(objects.symmetric_difference(other.objects)) != 0:
      return False
    if not np.isclose(self.value, other.value):
      if np.isnan(self.value) and np.isnan(other.value):
        pass
      else:
        return False
    return self.label == other.label

  def len(self):
    return len(self.objects)

  def union(self, other):
    """
    Computes the union of the current group with another group.
    :param Group other:
    """
    return Group(set(self.objects).union(other.objects))

  def get(self):
    return self.objects

  def copy(self):
    return Group(self.objects, value=self.value, label=self.label)

  @staticmethod
  def groupify(item):
    if isinstance(item, Group):
      return item
    else:
      try:
        return Group(item)
      except:
        raise ValueError("objects passed to Group is not iterable")


class GroupCollection(object):
  """Manage groups of similar objects."""

  def __init__(self, initial_groups=None):
    self.groups = set([])
    if initial_groups is not None:
      for group in initial_groups:
        self.add(Group.groupify(group))

  def __repr__(self):
    result = ""
    for group in self.groups:
      result = "%s\n%s" % (result, group)
    return result

  def add(self, group):
    self.groups.add(Group.groupify(group))

  def len(self):
    return len(self.groups)

  def intersection(self, other):
    """
    Constructs groups that are non-null intersections of groups.
    :param GroupCollection other:
    :return GroupCollection:
    Does not preserve value or label
    """
    new_collection = GroupCollection()
    for group in self.groups:
      for other_group in other.groups:
        new_group = group.intersection(other_group)
        if new_group.len() > 0:
          new_collection.add(new_group)
    return new_collection

  def unionDisjoint(self, other, **kwargs):
    """
    Constructs the union of disjoint groups, preserving group_values.
    :param GroupCollection other:
    :param dict kwargs: Arguments for GroupCollection constructor
    :return GroupCollection:
    """
    new_groups = set(self.groups).union(other.groups)
    if len(new_groups) < len(self.groups) + len(other.groups):
      raise ValueError("Not a disjoint union")
    return self.__class__(initial_groups=new_groups, **kwargs)

  def makeValueDF(self, default=np.nan, force_value=None):
    """
    :param object default: Value where none is specified
    :param float force_value: Force non-nan values to be this value
    :return pd.DataFrame: index is object, column is group label
        cells are the value for the group.
    """
    objects = self.flatten().get()
    num_rows = len(objects)
    num_cols = len(self.groups)
    # Construct a dataframe of nan's
    matrix = np.repeat(default, num_rows*num_cols)
    matrix = matrix.reshape(num_rows, num_cols)
    df = pd.DataFrame(matrix)
    # Assign values for groups
    columns = [g.id() for g in self.groups]
    columns.sort()
    df.columns = columns
    df.index = [str(o) for o in objects]
    for group in self.groups:
      if force_value is None:
        value = group.value
      else:
        value = force_value
      for obj in group.get():
        df.loc[str(obj), group.id()] = value
    return df

  def getGroupLabels(self):
    return [g.label for g in self.groups]

  def setGroupLabels(self, prefix=""):
    """
    Sets values of labels for all groups in the collection.
    """
    for idx, group in enumerate(self.groups):
      if len(prefix) > 0:
        stg = "%s-" % prefix
      else:
        stg = ""
      group.label = "%s%d" % (stg, idx)

  def equals(self, other):
    """
    Checks if this GroupCollection is the same as another.
    :param GroupCollection other:
    :return bool:
    """
    def isEqualSets(iter1, iter2):
      return len(set(set(iter1).symmetric_difference(iter2))) == 0
    #
    if not isEqualSets(self.getGroupLabels(), other.getGroupLabels()):
      return False
    if not isEqualSets(self.getGroupLabels(), other.getGroupLabels()):
      return False
    for group in self.groups:
      others = [g for g in other.groups if g.id() == group.id()]
      if len(others) != 1:
        return False
      other_group = others[0]
      if not np.isclose(group.value, other_group.value):
        return False
    return True

  def copy(self):
    # Construct the new groups
    return GroupCollection(
        initial_groups=[g.copy() for g in self.groups])

  def flatten(self):
    """
    :return list-str:
    """
    result = Group([])
    for group in self.groups:
      result = result.union(group)
    return result
    
