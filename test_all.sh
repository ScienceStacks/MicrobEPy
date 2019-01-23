#! /bin/bash
# Runs the all of the tests in the directory in python2 and python3
run_tests () {
  DIR=$1
  PATH=microbepy/tests/${DIR}
  for f in ${PATH}/test_*.py
    do
      echo ""
      echo "**** $f ****"
      $HOME/miniconda3/bin/python $f
    done
}


for d in common correlation data model plot search statistics
  do
    echo ""
    echo ">>>$d/<<<"
    echo ""
    run_tests $d
    echo ""
    echo ""
  done
