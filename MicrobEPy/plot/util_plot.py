"""Utilities for plotting."""

import constants as cn

import matplotlib.pyplot as plt

################################################
# Classes
################################################
class PlotParms(dict):
  """Container of plot parameters."""

  def __init__(self, is_initialize=True):
    """
    :param bool is_initialize: Set default values
    """
    super(self.__class__, self).__init__()
    if is_initialize:
      self.setattr(cn.PLT_FIGSIZE, (6, 4))
      self.setattr(cn.PLT_XLABEL, "Observed")
      self.setattr(cn.PLT_YLABEL, "Estimated")
      for attr in [cn.PLT_TITLE, cn.PLT_XLABEL, cn.PLT_YLABEL]:
        self.setattr(attr, '')

  def setattr(self, name, value):
    """
    Sets the value of the attribute if it is not present.
    """
    if not name in list(self.keys()):
      super(self.__class__, self).__setitem__(name, value)

  def isTrue(self, attr):
    """
    :param str attr:
    :return bool: True if attribute is present and True
    """
    if self.has_key(attr):
      return self.get(attr)
    else:
      return False

  def isFalse(self, attr):
    """
    :param str attr:
    :return bool: False if attribute is not present or False
    """
    return not self.isTrue(attr)

  def setTrueIfAbsent(self, attr):
    """
    :param str attr:
    """
    if not self.has_key(attr):
      self.setattr(attr, True)
    
  def do(self, is_plot=True):
    plt.xlabel(self[cn.PLT_XLABEL])
    plt.ylabel(self[cn.PLT_YLABEL])
    plt.title(self[cn.PLT_TITLE])
    # if no entry for PLT_LEGEND, take no action
    if self.has_key(cn.PLT_LEGEND):
      if self[cn.PLT_LEGEND] is None:
        plt.legend("")
      else:
        plt.legend(self[cn.PLT_LEGEND])
    if is_plot:
      plt.show()
      plt.close()
    if self.has_key(cn.PLT_XLIM):
      plt.xlim(self[cn.PLT_XLIM])
    if self.has_key(cn.PLT_YLIM):
      plt.ylim(self[cn.PLT_YLIM])