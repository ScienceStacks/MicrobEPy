#!/usr/bin/env python

from microbepy.common import constants as cn

import setuptools
import os
import shutil
import subprocess
import sys


DATA_DIR = "data_base"
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
    packages = [d for d in setuptools.find_packages() 
        if not ".tests" in d]
    setuptools.setup(name='microbepy',
        version='1.0',
        description='Python support for analysis of Microbial Communities',
        author='Joseph Hellerstein',
        author_email='jlheller@uw.edu',
        packages=packages,
        install_requires=getPipRequirements(),
        )
    # Other actions
    from microbepy.common import config
    config.setup(yaml_dict=cn.YAML_DEFAULT)
    print("--Conda installs")
    condaInstall()
  else:
    print("***No conda installation detected. Please install.")

if __name__ == '__main__':
  main()
