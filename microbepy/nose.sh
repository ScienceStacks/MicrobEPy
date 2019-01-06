# Run nosetests for all files in package
clear
nosetests -m pdb
nosetests --with-coverage --cover-inclusive --cover-package=microbepy
