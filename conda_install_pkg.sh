#!/bin/bash
# Runs a quiet install for a conda package
PKG=$1
#
conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda install $PKG
