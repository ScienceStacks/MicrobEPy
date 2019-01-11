#!/usr/bin/env python

import distutils.core
import os
import subprocess
import shutil


DATA_DIR = "data_base"
PROJECT_NAME = "microbepy"
DB_NAME = "%s.db" % PROJECT_NAME
REQUIREMENTS_FILE = "requirements.txt"
CONDA_PKGS = ['python', 'numpy', 'pandas', 'matplotlib',
'jupyter notebook', 'scikit-learn']

# Install required packages
def pipInstall():
  with open(REQUIREMENTS_FILE) as fd:
    requirements = fd.readlines()
  for requirement in requirements:
    pkg = requirement.replace('\n', '')
    subprocess.call(["pip", "install", pkg])

# Conda installs
def condaInstall():
  for pkg in CONDA_PKGS:
    subprocess.call(["conda", "install", pkg])

# Copy data to install directory
def copyData():
  curdir = os.path.realpath(os.curdir)
  dst = os.path.join(curdir, PROJECT_NAME)
  dst = os.path.join(dst, DATA_DIR)
  if not os.path.isdir(dst):
    os.mkdir(dst)
  dst = os.path.join(dst, DB_NAME)
  src = os.path.join(curdir, "Data")
  src = os.path.join(src, "data_model")
  src = os.path.join(src, DB_NAME)
  if not os.path.isfile(src):
    raise ValueError("Cannot file microbepy.db")
  shutil.copyfile(src, dst)

# Test to see if conda has been installed
def isCondaInstalled():
  result = True
  try:
    subprocess.call(["conda"], stdout=subprocess.DEVNULL)
  except:
    result = False
  return result

# Main logic
def main():
  if isCondaInstalled():
    copyData()
    distutils.core.setup(name='microbepy',
        version='1.0',
        description='Python support for analysis of Microbial Communities',
        author='Joseph Hellerstein',
        author_email='jlheller@uw.edu',
        packages=['microbepy'],
        package_dir={'microbepy': 'microbepy'},
        package_data={'microbepy': ['data_base/microbepy.db']},
        )
    print("--Pip installs.")
    pipInstall()
    print("--Conda installs.")
    condaInstall()
  else:
    print("***No conda installation detected. Please install.")

if __name__ == '__main__':
  main()
