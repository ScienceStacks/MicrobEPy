import os
import sys

PROJECT_DIRECTORY = "MicrobEPy"

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
    if path.find(PROJECT_DIRECTORY) < 0:
      path = last_path
      found = True
      break
  if not found:
    raise RuntimeError("Could not find project path.")
  return path

path = os.path.join(getProjectDirectory(), PROJECT_DIRECTORY)
sys.path.append(path)  # Get the directory for python code
from project_base import addPythonPaths
addPythonPaths()
