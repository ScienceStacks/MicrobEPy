# Run nosetests for all files in package
clear
nosetests -m pdb
nosetests --with-coverage --cover-inclusive --cover-package=MicrobEPy --cover-erase -m pdb
