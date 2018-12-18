"""Calculate Equivalence Classes."""

class EquivalenceClass(object):

  def __init__(self, iterable, relation):
    """
    :param iterable iterable: collection of objects to be partitioned
    :param BooleanBinaryFunction relation:
        relation: equivalence relation. I.e. relation(o1,o2) evaluates to True
            if and only if o1 and o2 are equivalent
    """
    self._iterable = iterable
    self._relation = relation
    self.classes = []  # Collection of equivalent items
    self.partitions = {}  # Relates item to its set
    self.ids = []  # List of identifiers associated witht the corresponding class

  def do(self):
    """
    Construct the classes, partitions and ids.
    """
    self.classes, self.partitions, self.ids = self._equivalence_enumeration()

  def _equivalence_partition(self):
    """Partitions a set of objects into equivalence classes
    Returns: classes, partitions
        classes: A sequence of sets. Each one is an equivalence class
        partitions: A dictionary mapping objects to equivalence classes
    """
    classes = []
    partitions = {}
    for o in self._iterable:  # for each object
      # find the class it is in
      found = False
      for c in classes:
        if self._relation(next(iter(c)), o):  # is it equivalent to this class?
          c.add(o)
          partitions[o] = c
          found = True
          break
      if not found:  # it is in a new class
        classes.append(set([o]))
        partitions[o] = classes[-1]
    return classes, partitions

  def _equivalence_enumeration(self):
    """
    Returns: classes, partitions, ids
        classes: A sequence of sets. Each one is an equivalence class
        partitions: A dictionary mapping objects to equivalence classes
        ids: A dictionary mapping objects to the indices of their equivalence classes
    """
    classes, partitions = self._equivalence_partition()
    ids = {}
    for i, c in enumerate(classes):
      for o in c:
        ids[o] = i
    return classes, partitions, ids

  def validate(self):
    """
    Validates the partitions obtained.
    """
    for o, c in self.partitions.items():
        for _c in self.classes:
          predicate = (o in _c) ^ (not _c is c)
          if not predicate:
            raise ValueError("Invalid partitions")
