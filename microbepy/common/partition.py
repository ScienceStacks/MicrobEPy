"""Manipulate partitions of a set of elements."""

import collections
import copy

SetMutation = collections.namedtuple('SetMutation',
    'cur_src new_src cur_dst new_dst')


class Partition(object):

  def __init__(self, elements, partition=None):
    """
    :param enumerable-of-object elements: Must be unique
    :param set of sets: partition of elements. If None,
        initialized to singletons
    """
    self._validateSet(elements)
    self.elements = set(copy.deepcopy(elements))
    if partition is None:
      partition = [ set(e) for e in elements]
    else:
      self._validateSet(partition)
      items = []
      [items.extend(s) for s in partition]
      if not set(items) == self.elements:
        raise ValueError("Invalid partition specified.")
      partition = copy.deepcopy(partition)
    self.sets = partition

  def _validateSet(self, a_set):
    for ele in a_set:
      if list(a_set).count(ele) > 1:
        raise ValueError("Not a set.")
    
  def isPresent(self, a_set):
    """
    :return bool: True if a_set is present
    """
    return a_set in self.sets

  def _findSet(self, a_set):
    """
    :param list-of-sets sets:
    :param set a_set:
    :return int: index of the set or -1
    """
    try:
      idx = self.sets.index(a_set)
    except ValueError:
      idx = -1
    return idx
  
  def move(self, element, cur_dst, is_update=False):
    """
    Moves the element to the destination set.
    If the partition is to be updated (is_update),
    (a) if the source set is empty, it is deleted;
    (b) if the destination set is not in the partition,
    it is added.
    :param object element:
    :param set dst:
    :param bool is_update: Update the sets if True
    :return SetMutation:
    """
    cur_src = self.findSetWithElement(element)
    new_src = cur_src.difference([element])
    new_dst = cur_dst.union([element])
    result = SetMutation(cur_src=cur_src, new_src=new_src,
        cur_dst=cur_dst, new_dst=new_dst)
    if is_update:
      # Source set
      if len(new_src) == 0:
        self.sets.remove(cur_src)
      else:
        idx = self._findSet(cur_src)
        self.sets[idx] = new_src
      # Destination set
      if len(cur_dst) == 0:
        self.sets.append(new_dst)
      else:
        idx = self._findSet(cur_dst)
        self.sets[idx] = new_dst
    return result

  def findSetWithElement(self, element, 
      is_present=True):
    """
    Finds the set with specified element.
    :param object element:
    :param bool is_present: if True, must be present
    :return set:
    """
    sets = [s for s in self.sets if element in s]
    if len(sets) > 1:
      raise RuntimeError("Should have at most instance of values.")
    if len(sets) == 0:
      if is_present:
        raise ValueError("%s is not present" % str(element))
      result = cn.NULL_SET 
    else:
      result = sets[0]
    return result
