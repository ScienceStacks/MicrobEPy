#!/bin/bash
# Set-up miniconda
PY=3
CONDA=$HOME/miniconda${PY}/bin/conda
rm -rf $HOME/miniconda${PY}
wget https://repo.continuum.io/miniconda/Miniconda${PY}-latest-Linux-x86_64.sh
bash Miniconda${PY}-latest-Linux-x86_64.sh
rm Miniconda${PY}-latest-Linux-x86_64.sh
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
