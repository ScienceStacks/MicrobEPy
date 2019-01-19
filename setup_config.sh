#!/bin/bash
mkdir ../.microbepy
cd Data/data_model
PATH="`pwd`/microbepy.db"
cd ../..
echo "SQLDB_PATH: ${PATH}" > ../.microbepy/config.yml
/bin/cat ../.microbepy/config.yml
