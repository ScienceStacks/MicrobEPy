import os
import sys
import warnings

# Get rid of bogus warning messages from sklearn
def warn(*args, **kwargs):
  pass
warnings.warn = warn


PROJECT_NAME = "microbepy"

def getProjectDirectory():
  """
  Finds the path to the folder containing this project.
  :return str path:
  """
  curdir = os.getcwd()
  paths = []
  # Find the list of subpaths to this name
  while len(curdir) > 1:
    paths.append(curdir)
    curdir = os.path.split(curdir)[0]
  paths.reverse()
  # Find the path to the directory for this project
  found = False
  for path in paths:
    if PROJECT_NAME in path:
      found = True
      break
  if not found:
    raise RuntimeError("Could not find project path.")
  return path

project_dir = getProjectDirectory()
main_code_path = os.path.join(project_dir, PROJECT_NAME)
sys.path.append(main_code_path)
