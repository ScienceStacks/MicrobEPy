"""Iteratively supplies combinations with optional exclusions."""

import __init__
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
    if excludes is None:
      excludes = []
    self._excludes = excludes
    self._combination_iterator = None
    self._current_size = 0

  def __iter__(self):
    while True:
      if self._combination_iterator is None:
        self._current_size += 1
        if self._current_size > self._max_size:
          break
        self._combination_iterator = itertools.combinations(
            self._elements, self._current_size)
      for vals in self._combination_iterator:
        vals = [v for v in vals]
        if not vals in self._excludes:
          yield vals
      self._combination_iterator = None
