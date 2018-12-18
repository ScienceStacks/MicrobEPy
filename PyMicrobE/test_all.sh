#! /bin/bash
# Runs the all of the tests in the directory in python2 and python3
run_tests () {
  for f in test_*.py
    do
      echo ""
      echo "**** $f ****"
      python $f
      #python3 $f
    done
}


for d in common correlation data data_access model plot search statistics
  do
    cd $d
    echo ""
    echo ">>>$d/<<<"
    echo ""
    run_tests
    cd ..
    echo ""
    echo ""
  done
