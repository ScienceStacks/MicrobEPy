"""Iteratively supplies combinations with optional exclusions."""

import __init__
import util
import constants as cn

import itertools


class CombinationIterator(object):

  def __init__(self, elements, max_size, excludes=None):
    """
    :param list-object elements: Elements to construct combinations
    :param int max_size: Maximum number of elements in a combination
    :param list-object excludes: Combinations to be excluded
    """
    self._elements = elements
    self._max_size = max_size
    self._excludes = util.setNoneList(excludes)
    self._combination_iterator = None
    self._current_size = 0

  def __iter__(self):
    return self

  def next(self):
    while True:
      if self._combination_iterator is None:
        self._current_size += 1
        if self._current_size > self._max_size:
          raise StopIteration
        self._combination_iterator = itertools.combinations(
            self._elements, self._current_size)
      for vals in self._combination_iterator:
        vals = [v for v in vals]
        if not vals in self._excludes:
          return vals
      self._combination_iterator = None
