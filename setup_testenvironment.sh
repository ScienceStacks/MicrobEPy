# Creates a conda environment in which to run the tests
conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda info -a

if [ ! -e /home/ubuntu/miniconda3/envs/test-environment ]; then
  echo "***Create environment***"
  conda create -q -n test-environment python=3.6 numpy pandas matplotlib jupyter scikit-learn
  pip install xlrd
  pip install nose
  pip install coverage
fi
conda info --envs
echo "***Done!."
echo "***Next: 'source activate test-environment'"
echo "***Then: cd MicrobeEPy; nosetests"
