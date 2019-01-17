"""
Manages configuration information.
Configuration information is maintained in the .microbepy
directory of the account root in the file config.yaml
"""

from microbepy.common import constants as cn

import os
import yaml

YAML_DEFAULT = {cn.SQLDB_PATH_NAME: cn.SQLDB_PATH}

def initialize():
  """
  Sets up the directory if needed
  """
  if os.path.isdir(cn.CONFIG_DIR_PATH):
    return
  os.mkdir(cn.CONFIG_DIR_PATH)

def setup(yaml_default=YAML_DEFAULT, is_forced=False):
  """
  Ensures that the cn.CONFIG_FILE is present.
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
  with open(cn.CONFIG_FILE_PATH, 'w') as fd:
    fd.writelines(lines)
  return yaml_dict
