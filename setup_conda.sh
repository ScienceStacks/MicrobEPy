#!/bin/bash
# Set-up miniconda
PY=3
CONDA=$HOME/miniconda${PY}/bin/conda
rm -rf $HOME/miniconda${PY}
wget https://repo.continuum.io/miniconda/Miniconda${PY}-latest-Linux-x86_64.sh
bash Miniconda${PY}-latest-Linux-x86_64.sh
rm Miniconda${PY}-latest-Linux-x86_64.sh
