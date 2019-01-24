"""
Manages configuration information.
Configuration information is maintained in the .microbepy
directory of the account root in the file cn.CONFIG_FILE
"""

from microbepy.common import constants as cn

import os
import yaml

def initialize():
  """
  Sets up the directory if needed
  """
  if os.path.isdir(cn.CONFIG_DIR_PATH):
    return
  os.mkdir(cn.CONFIG_DIR_PATH)

def setup(yaml_dict=cn.YAML_DEFAULT):
  """
  Sets up the configuration file if it doesn't exist.
  :param str yaml_dict: yaml used if none is present
  """
  initialize()
  if os.path.isfile(cn.CONFIG_FILE_PATH):
    pass
  else:
    with open(cn.CONFIG_FILE_PATH, 'w') as fd:
      yaml.dump(yaml_dict, fd, default_flow_style=False)

def get(key=None):
  """
  Returns the value of a configuration key.
  :param str key: Configuration key. If key
    is None, then returns entire dictionary.
  :return value:
  :Raises KeyError: Key not present
  """
  if os.path.isfile(cn.CONFIG_FILE_PATH):
    with open(cn.CONFIG_FILE_PATH, 'r') as fd:
      yaml_dict = yaml.load(fd)
      if key is None:
        result = yaml_dict
      else:
        result = yaml_dict[key]
  else:
    raise KeyError("%s does not exist." % key)
  return result
