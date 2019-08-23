"""Utilities for plotting."""

from microbepy.common import constants as cn

import matplotlib.pyplot as plt

################################################
# Classes
################################################
class PlotParms(dict):
  """Container of plot parameters."""

  def __init__(self, is_initialize=True,
      fontsize_label=12, fontsize_title=14):
    """
    :param bool is_initialize: Set default values
    """
    super(self.__class__, self).__init__()
    self.fontsize_label = fontsize_label
    self.fontsize_title = fontsize_title
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
    if attr in self:
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
    if not attr in self:
      self.setattr(attr, True)
    
  def do(self, is_plot=True):
    plt.xlabel(self[cn.PLT_XLABEL],
        fontsize=self.fontsize_label)
    plt.ylabel(self[cn.PLT_YLABEL],
        fontsize=self.fontsize_label)
    plt.title(self[cn.PLT_TITLE], fontsize=self.fontsize_title)
    # if no entry for PLT_LEGEND, take no action
    if cn.PLT_LEGEND in self:
      if self[cn.PLT_LEGEND] is None:
        plt.legend("")
      else:
        plt.legend(self[cn.PLT_LEGEND])
    if is_plot:
      plt.show()
      plt.close()
    if cn.PLT_XLIM in self:
      plt.xlim(self[cn.PLT_XLIM])
    if cn.PLT_YLIM in self:
      plt.ylim(self[cn.PLT_YLIM])
