#!/bin/bash
# Cleans up files used in testing
DIRS="data_access phenotype_explorer regression statistics growth_data ."

mv_file () {
  FILE=$1
  if [ -e $FILE ]
    then mv $FILE /tmp
  fi
}

clean_dir() {
  cd $1
  mv_file api.pcl
  mv_file *.pyc
  mv_file *_test.csv
  cd ..
}

for d in $DIRS
  do
    clean_dir $d
  done
