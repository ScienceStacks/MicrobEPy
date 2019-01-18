#!/usr/bin/env python

from microbepy.common import config
from microbepy.common import constants as cn

from setuptools import setup, find_packages
import os
import shutil
import subprocess
import sys


DATA_DIR = "data_base"
PROJECT_NAME = "microbepy"
DB_NAME = "%s.db" % PROJECT_NAME
REQUIREMENTS_FILE = "requirements.txt"
CONDA_PKGS = ['python', 'numpy', 'pandas', 'matplotlib',
'jupyter notebook', 'scikit-learn']

# Install required packages
def getPipRequirements():
  with open(REQUIREMENTS_FILE) as fd:
    requirements = fd.readlines()
  return requirements

# Conda installs
def condaInstall():
  for pkg in CONDA_PKGS:
    subprocess.call(["conda_pkg_install.sh"])

# Test to see if conda has been installed
def isCondaInstalled():
  result = True
  try:
    subprocess.call(["conda"], stdout=subprocess.DEVNULL)
  except:
    result = False
  return result

# Conda installs
def condaInstall():
  for pkg in CONDA_PKGS:
    subprocess.call(["conda", "install", pkg])

# Main logic
def main():
  if isCondaInstalled():
    # Setup the configuration directory
    sqldb_path = os.getcwd()
    for ele in ["Data", "data_model", cn.SQLDB_FILE]:
      sqldb_path = os.path.join(sqldb_path, ele)
    yaml_settings = {cn.SQLDB_PATH_NAME: sqldb_path}
    config.setup(yaml_default=yaml_settings, is_forced=True)
    packages = [d for d in find_packages() if not ".tests" in d]
    setup(name='microbepy',
        version='1.0',
        description='Python support for analysis of Microbial Communities',
        author='Joseph Hellerstein',
        author_email='jlheller@uw.edu',
        packages=packages,
        install_requires=getPipRequirements(),
        package_data={'microbepy': ['data_base/microbepy.db']},
        )
    print("--Conda installs")
    condaInstall()
  else:
    print("***No conda installation detected. Please install.")

if __name__ == '__main__':
  main()
