"""
Manages configuration information.
Configuration information is maintained in the .microbepy
directory of the account root in the file config.yaml
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

def setup(yaml_default=cn.YAML_DEFAULT, is_forced=False):
  """
  Ensures that the cn.SQLDB_FILE is present.
  :param str yaml_default: yaml used if none is present
  :param bool is_forced: force the yaml file to the default
  """
  initialize()
  if os.path.isfile(cn.CONFIG_FILE_PATH) and (not is_forced):
    with open(cn.CONFIG_FILE_PATH, 'r') as fd:
      yaml_dict = yaml.load(fd)
  else:
    yaml_dict = yaml_default
  lines = yaml.dump(yaml_dict, default_flow_style=False)
  if os.path.isdir(cn.CONFIG_DIR_PATH):
    with open(cn.CONFIG_FILE_PATH, 'w') as fd:
      yaml.dump(yaml_dict, fd, default_flow_style=False)
  return yaml_dict
