"""Codes that establish the base directory structure."""

import os
import sys


PROJECT_ROOT = "MicrobEPy"
PYTHON_SUBDIRECTORIES = [
    "statistics", "model", "data_access", "correlation",
    "data", "plot", "search", "common",
    ]


def addPythonPaths():
  """
  Adds the paths needed for python code.
  """
  project_dir = getProjectDirectory()
  # Directory of python codes
  main_code_path = os.path.join(project_dir, PROJECT_ROOT)
  # Directory of python codes
  sys.path.append(main_code_path)
  for directory in PYTHON_SUBDIRECTORIES:
    path = os.path.join(main_code_path, directory)
    sys.path.append(path)

def getProjectDirectory():
  """
  :return str path:
  """
  path = os.getcwd()
  # Go up to coevolution
  max_iteration = path.count("/")
  found = False
  for n in range(max_iteration):
    last_path = path
    path = os.path.split(path)[0]
    if path.find(PROJECT_ROOT) < 0:
      path = last_path
      found = True
      break
  if not found:
    raise RuntimeError("Could not find project path.")
  return path
