# Configuration constants for MicrobEPy

import os

# Path to the SQL data.
# If None and in a github repo, uses Data/data_model/microbepy.db
possible_db = os.path.realpath('data_base')
possible_db = os.path.join(possible_db, 'microbepy.db')
if os.path.isfile(possible_db):
  SQLDB_PATH = possible_db
else:
  SQLDB_PATH = None
