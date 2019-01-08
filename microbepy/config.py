"""Functions to update configurations"""

import yaml

CONFIG_FILE = "config.yaml"

def get(key=None, config_file=CONFIG_FILE):
  """
  Returns value for key. Otherwise, returns entire data map.
  """
  with open(config_file) as f:
    data_map = yaml.safe_load(f)
  data_map = {k: None if v == 'None' else v 
      for k, v in data_map.items()}
  if key is not None:
    value = data_map[key]
  else:
      value = data_map
  return value

def set(key, value, config_file=CONFIG_FILE):
  """
  Sets the value of a configuration parameter.
  """
  try:
    data_map = get(config_file=config_file)
  except FileNotFoundError:
    data_map = {}
  import pdb; pdb.set_trace()
  data_map = {k: 'None' if v is None else v 
      for k, v in data_map.items()}
  data_map[key] = value
  with open(config_file, "w") as f:
    yaml.dump(data_map, f)
