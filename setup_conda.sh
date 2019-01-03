#!/bin/bash
# Does setup after miniconda is installed
CONDA=$HOME/miniconda3/bin/conda
# Do Conda installs
$CONDA update -n base conda
$CONDA install python=3.6.4  # Latest release for Tellurium
$CONDA install numpy
$CONDA install pandas
$CONDA install matplotlib
$CONDA install jupyter notebook
$CONDA install scikit-learn
#
sudo apt install python-setuptools
# Other pip installs
pip install xlrd
pip install nose
pip install coverage
