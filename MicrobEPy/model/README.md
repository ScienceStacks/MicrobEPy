# Directory Overview

This directory addresses model building. Two types of models
are considered: regression and classification. Models are
evaluated using cross validation, a process by which a model
is trained on one subset of the data and then scored using
the remaining data. 

Two kinds of minds are considered. Isolate models predict
phenotype from isolate characteristics such as phenotype
from ancestral pairings. Genome models predict phenotype from
mutations.

- cv_\* are cross valdiation codes
- isolate_\* files are isolate model codes
- *_model* files are abstract classes for regression and
  classification models
- test_\* files are test files
- \*_transformer files transform data that are input to files,
  either predictors or dependent variable
