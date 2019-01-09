#!/usr/bin/env python

from distutils.core import setup
import os
import shutil


DATA_DIR = "data_base"
PROJECT_NAME = "microbepy"
DB_NAME = "%s.db" % PROJECT_NAME

# Copy data to install directory
def copyData():
  curdir = os.path.realpath(os.curdir)
  dst = os.path.join(curdir, PROJECT_NAME)
  dst = os.path.join(dst, DATA_DIR)
  if not os.path.isfile(dst):
    os.mkdir(dst)
  dst = os.path.join(dst, DB_NAME)
  src = os.path.join(curdir, "Data")
  src = os.path.join(src, "data_model")
  src = os.path.join(src, DB_NAME)
  if not os.path.isfile(src):
    raise ValueError("Cannot file microbepy.db")
  shutil.copyfile(src, dst)

if __name__ == '__main__':
  copyData()
  setup(name='microbepy',
      version='1.0',
      description='Python support for analysis of Microbial Communities',
      author='Joseph Hellerstein',
      author_email='jlheller@uw.edu',
      packages=['microbepy'],
      )
