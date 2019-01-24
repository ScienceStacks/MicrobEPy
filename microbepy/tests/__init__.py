# Eliminate warnings from tests of all modules

import warnings

def dummy_warn(*pargs, **kwargs):
  pass

warnings.warn = dummy_warn
