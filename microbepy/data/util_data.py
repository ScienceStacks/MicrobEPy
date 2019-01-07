"""Utilities for creating data."""

import microbepy_init
import constants as cn
from isolate_regression import IsolateRegression
from isolate import Isolate

import numpy as np
import pandas as pd
import scipy.stats as stats
import util

MAX_STD = 3.0
MUTATION_GROUP_STRING = "--"  # Mutation group string

def makeCultureIsolateMutationDF(is_separate_species=True):
  """
  Retrieves the data required for making plots
  :return pd.DataFrame:
    cn.LINE, cn.KEY_CULTURE, KEY_ISOLATE_DVH, KEY_ISOLATE_MMP
    cn.KEY_MUTATION, cn.GGENE_ID, cn.GENE_ID, cn.POSITION
    cn.RATE, cn.YIELD
  """
  query1 = '''
  select distinct key_culture, line, key_isolate, 
      gene_id, ggene_id, position, key_mutation, effect,
      rate, yield
    from genotype_phenotype
      where species_mix = 'B' 
          and is_an_mutation = 0
  '''
  df1 = util.readSQL(query1)
  if is_separate_species:
    # Separate columns for species
    del df1[cn.KEY_ISOLATE]
    query2 = '''
    select distinct key_culture,
            key_isolate as key_isolate_dvh, 
            key_isolate_mmp 
          from genotype_phenotype,
            (select distinct key_isolate as key_isolate_mmp, 
                line as line_mmp,
                key_culture as key_culture_mmp from genotype_phenotype 
              where species='M' 
                and line_mmp != 'AN') sub 
          where species_mix = 'B' 
              and species = 'D' 
              and key_culture_mmp = key_culture
              and line != 'AN'
    '''
    df2 = util.readSQL(query2)
    sel = [Isolate.create(i).experiment == cn.EXPERIMENT_CI 
          for i in df2[cn.KEY_ISOLATE_DVH]]
    df2 = df2.loc[sel]
    df = df1.merge(df2, on=cn.KEY_CULTURE, how='inner')
  else:
    sel = [Isolate.create(i).experiment == cn.EXPERIMENT_CI 
          for i in df1[cn.KEY_ISOLATE]]
    df1 = df1.loc[sel]
    df = df1
  df[cn.POSITION] = df[cn.POSITION].apply(lambda v: int(v))
  return df

def filterOutlierCultures(df):
  """
  Removes cultures that have large residuals for the
  NonParamatericIsolateRegression.
  :param pd.DataFrame df:  has column cn.KEY_CULTURE
  :return pd.DataFrame:
  """
  cultures = IsolateRegression.getSmallResidualCultures(
      max_std=MAX_STD)
  #
  rows = []
  for _, row in df.iterrows():
    if row[cn.KEY_CULTURE] in cultures:
      rows.append(row)
  return pd.DataFrame(rows)

def makeGausianMixture(means, stds, num_rows, num_extracols,
    prob=0.5):
  """
  Creates predictor and dependent variable dataframes. 
  There is a predictor column for each mean and each 
  extra column. Predictor values are binary and are assigned
  randomly. The
  n-th dependent variable is a gausian mixture of the means
  that have a binary 1 in their predictor column.
  :param list-float means:
  :param list-float stds:
  :param int num_rows: number of rows in result
  :param int num_extracols: number of extra columns that do
                            not affect the dependent variable
  :param float prob: probability of a 1 in df_X
  :return pd.DataFrame, pd.DataFrame: df_X, df_y
  """
  COL_Y = 'y'
  num_cols = len(means) + num_extracols
  # construct the predictors
  df_X = pd.DataFrame(
       stats.bernoulli.rvs(prob, size=(num_rows, num_cols)))
  predictor_columns = df_X.columns.tolist()[0:len(means)]
  extra_columns = df_X.columns.tolist()[len(means):]
  for idx, col in enumerate(predictor_columns):
    df_X.rename(columns={col: "X_%2.4f" % means[idx]}, inplace=True)
  for idx, col in enumerate(extra_columns):
    df_X.rename(columns={col: "E_%d" % idx}, inplace=True)
  predictor_columns = df_X.columns.tolist()[0:len(means)]
  extra_columns = df_X.columns.tolist()[len(means):]
  # Construct the dependent variables
  df_y = pd.DataFrame({
    COL_Y: np.repeat(0, num_rows)
    })
  for idx, col in enumerate(predictor_columns):
    sel = df_X[col] == 1
    values = np.random.normal(means[idx], stds[idx], num_rows)
    y_values = [y + v if s else y 
             for s, v, y in zip(sel, values, df_y[COL_Y])]
    df_y[COL_Y] = y_values
  #
  return df_X, df_y

def makeNoisyOr(predictor_probs, num_rows, noise_prob):
  """
  Creates a dependent variable (y) that is the OR of
  predictor variables (x_i) with noise added. The resulting
  binary variable has the value of the OR with probability
  1 - p_noise and the value of a Bernoulli(0.5) with probability
  p_noise. The dependent variables is distributed as a
  Bernoulli(0.5). P(y | x_i, p_i) = 1, where p_i is the fraction
  of the data in which y_i = x_i.
  :param list-float predictor_probs: probability distribution
                           for the predictor variables.
  :param int num_rows: number of rows in result
  :param float noise_prob:
  :return pd.DataFrame, pd.DataFrame, float: df_X, df_y, score
    score - maximum accuracy achievable by a classifier
  """
  def concat(mat1, mat2):
    if mat1 is None:
      return mat2
    if mat2 is None:
      return mat1
    return np.concatenate((mat1, mat2), axis=1)
  def randomizeIndex(dfs):
    """
    Randomizes a set of dataframes in the same way.
    :param list-DataFrame dfs:
    :return list-DataFrame:
    """
    indices = np.random.permutation(range(len(dfs[0])))
    results = []
    for df in dfs:
      df[cn.INDEX] = indices
      df = df.sort_values(cn.INDEX)
      del df[cn.INDEX]
      df.index = range(len(df))
      results.append(df)
    return results
  #
  # Check the inputs
  if not np.isclose(np.sum(predictor_probs), 1.0):
    msg = "Predictor probabilities must be a distribution."
    raise ValueError(msg)
  # Constants used
  col_y = 'y'
  num_cols = len(predictor_probs)
  rand_prob = 0.5
  num_nonnoise_rows = int(np.round(num_rows*(1 -noise_prob)))
  num_noise_rows = num_rows - num_nonnoise_rows
  # Construct the dependent variable
  ys = stats.bernoulli.rvs(rand_prob, size=(num_rows, 1))
  # construct the predictors
  columns = ["P_%d" % n for n,_ in enumerate(predictor_probs)]
  dfs = []
  cur_row = 0
  for idx, prob in enumerate(predictor_probs):
    cur_num = int(np.round(num_rows*prob*(1 -noise_prob)))
    # First block
    if idx > 0: 
      mat1 = stats.bernoulli.rvs(rand_prob, size=(cur_num, 
         idx))
      mat1.reshape(idx, cur_num)
    else:
      mat1 = None
    # Predicator variable equal to y
    mat2 = np.array(ys[cur_row:(cur_row+cur_num)])
    mat2.reshape(cur_num, 1)
    # Third block
    if num_cols - idx > 1:
      num_cols3 = num_cols - idx - 1
      mat3 = stats.bernoulli.rvs(rand_prob, size=(cur_num, num_cols3))
      mat3.reshape(num_cols3, cur_num)
    else:
      mat3 = None
    # Assemble this section
    mat = concat(mat1, mat2)
    mat = concat(mat, mat3)
    df = pd.DataFrame(mat, columns=columns)
    dfs.append(df)
    cur_row += cur_num
  # Add the remaining rows
  rand_mat = stats.bernoulli.rvs(rand_prob, 
      size=(num_noise_rows, num_cols))
  df = pd.DataFrame(rand_mat, columns=columns)
  dfs.append(df)
  # Calculate the maximum accuracy of a classifier of these data
  # The accuracy is the rand_prob for the rows in which predictor
  # values are assigned randomly. The accuracy is 1 where one
  # predictor has a value equal to the dependent variable.
  score = rand_prob*noise_prob + 1 - noise_prob
  #
  df_X = pd.concat(dfs, sort=True)
  df_y = pd.DataFrame(ys, columns=[col_y])
  (df_X, df_y)  = randomizeIndex((df_X, df_y))
  #
  return df_X, df_y, score

def makeNoisyBin(num_rows, noise_prob, 
    predicate=lambda v: stats.bernoulli.rvs(0.5, size=1)[0]):
  """
  Predictor values are digits of a binary representation of the
  index. The dependent variable is in {0, 1}.
  A fraction noise_prob of dependent variables are randomly selected
  their predictor variables are duplicated and the dependent
  variable is negated. This is done in a way so that
  Half of the y values are 1s.
  :param int num_rows: number of rows in result
  :param float noise_prob:
  :param FunctionOfInt predicate: returns 0 or 1
  :return pd.DataFrame, pd.DataFrame: 
    df_X columns: P_*
    df_y column: cn.VALUE
  The maximum accuracy of a classifier is: 1/(1 + f), where f
  is the noise fraction.
  """
  def selectRandomIndexForValue(df_y, val, count):
    """
    Flips the value of randomly selected dependent variables
    """
    df = df_y[df_y[cn.VALUE] == val]
    return np.random.permutation(df.index.tolist())[0:count]
  #
  indices = range(num_rows)
  num_columns = int(np.ceil(np.log2(num_rows)))
  # Create the initial dependent variable
  df_y = pd.DataFrame({cn.VALUE: indices}, index=indices)
  df_y = df_y.applymap(lambda v: predicate(v))
  ys = df_y[cn.VALUE]
  frac_ones = (1.0*sum(ys))/len(ys)
  # Create the predictor variables
  xs = [ [int(v) for v in util.extendBin(n, num_columns)]
        for n in indices]
  columns = [("P_%d" % n) for n in range(num_columns)]
  columns.reverse()  # use conventional order of binary digits
  df_X = pd.DataFrame(xs, columns=columns, index=indices)
  # Add noise as required.
  num_random = np.round(num_rows*noise_prob)
  num_random_dict = {
      cn.ONE: int(np.round(frac_ones*num_random)),
      cn.ZERO: int(np.round((1.0-frac_ones)*num_random)),
      }
  if num_random_dict[cn.ONE] > 0:
    for val in cn.BINARY_VALUES:
      indices = selectRandomIndexForValue(df_y, val,
          num_random_dict[val])
      df_y = df_y.append(df_y.loc[df_y.index[indices]])
      df_y.index = range(len(df_y))
      df_X = df_X.append(df_X.loc[df_X.index[indices]])
      df_X.index = range(len(df_X))
  #
  return df_X, df_y

def makeGroups(num_rows, columns=None):
  """
  Creates a two column grouping dataframe used by splitter.
  :param int num_rows:
  :param list-str columns: columns in grouping dataframe
  :return pd.DataFrame:
  """
  def createValues(df, name):
    df[name] = [name + str(n) for n in range(num_rows)]
    return df
  #
  if columns is None:
    columns = ['a', 'b']
  df = pd.DataFrame()
  for col in columns:
    createValues(df, col)
  return df

def makeClassificationData(df_X, df_y, column, percentile):
  """
  Transforms a continuous dependent variable in df_y into a
  categorical variable based on a percentile threshold.
  Classes are 0, 1 based on the percentile criteria.
  Adjusts df_X where values are dropped.
  :param pd.DataFrame df_X:
  :param pd.DataFrame df_y:
  :param str column: column in df_y
  :param float percentile: percentile used to classify data
  """
  def makeDropSelectors():
    """
    Selectors are indices.
    :return list-int, list-int, list-int:
      sel0 - selected for class 0
      sel1 - selected for class 1
      sel_drop - indices to drop
    """
    value0 = np.percentile(df_y[column], percentile)
    sel0 = [n for n,v in enumerate(df_y[column]) if v <= value0]
    value1 = np.percentile(df_y[column], 100 - percentile)
    sel1 = [n for n,v in enumerate(df_y[column]) if v > value1]
    sel_drop = set(range(len(df_y))).difference(sel0)
    sel_drop = set(sel_drop).difference(sel1)
    return sel0, sel1, sel_drop
  #
  sel0, sel1, sel_drop = makeDropSelectors()
  ser = df_y[column]
  ser.iloc[sel0] = 0
  ser.iloc[sel1] = 1
  df_y[column] = ser
  y_drops = [df_y.index[n] for n in sel_drop]
  X_drops = [df_X.index[n] for n in sel_drop]
  [df_y.drop(d, inplace=True) for d in y_drops]
  [df_X.drop(d, inplace=True) for d in X_drops]

def makeMutationGroupDF(df_X, df_y, mutations):
  """
  Assigns isolates to groups based on their mutations.
  A group is a binary encoding of the presenence (binary 1)
  or absence (binary 0) of a mutation. The group is the
  decimal value of this binary number.
  :param pd.DataFrame df_X: 
      indexed by cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP
      columns of mutations
  :return pd.DataFrame: 
      cn.KEY_ISOLATE_DVH, cn.KEY_ISOLATE_MMP, cn.LINE
      cn.GROUP - int whose binary value ecodes the
                 presence/absence of mutations
      cn.DEPVAR - value of the dependent variable
  """
  df_result = df_X[mutations].copy()
  util.addBinColumn(df_result)  #  Add the binary valued column
  df_result.rename(columns={cn.VALUE: cn.GROUP}, inplace=True)
  for mutation in mutations:
    del df_result[mutation]
  df_result = df_result.reset_index()
  df_result[cn.LINE] = [Isolate.create(r[cn.KEY_ISOLATE_DVH]).line
      for _,r in df_result.iterrows()]
  col_y = df_y.columns.tolist()[0]
  df_result[cn.DEPVAR] = df_y[col_y].tolist()
  return df_result

def generatePhenotypeDataForMutations(df_X, df_y, mutations, num_replications):
  """
  Generates instances of data in groups specified by the mutations.
  :param pd.DataFrame df_X: mutation matrix
  :param pd.DataFrame df_y: phenotype data
  :param list-str mutations: list of mutations considered
  :param num_replications: number of replications of data instances
  :return pd.DataFrame:
     cn.REPLICATION, cn.GROUP, cn.VALUE
  """
  df_group = makeMutationGroupDF(df_X, df_y, mutations)
  groups = util.getFirstColumn(df_group.groupby(cn.GROUP).count())
  #
  return generatePhenotypeData(util.getFirstColumn(df_y), 
      groups, num_replications)

def generatePhenotypeData(values, group_sizes, num_replications):
  """
  Generates instances of data in groups specified by the mutations.
  :param list-float values: phenotype data
  :param list-int group_sizes: count of values to assign in a group
  :param num_replications: number of replications of data instances
  :return pd.DataFrame:
     cn.REPLICATION, cn.GROUP, cn.VALUE
  :raises ValueError:
  """
  if len(values) != sum(group_sizes):
    raise ValueError("Invalid specification of group_sizes. Must sum to len(values)")
  #
  # Construct basic dataframe
  df_result = pd.DataFrame({cn.VALUE: values})
  group_vals = []
  [[group_vals.append(m) for _ in range(s)]
      for m, s in enumerate(group_sizes)]
  df_result[cn.GROUP] = group_vals
  df_result = pd.concat([df_result]*num_replications,
      ignore_index=True, sort=True)
  # Construct the replications
  replications = []
  [replications.extend(np.repeat(n, len(values)))
      for n in range(num_replications)]
  df_result[cn.REPLICATION] = replications
  # Randomly order the y values within each replication,
  # keeping the same number of values in a group.
  groups = df_result[cn.GROUP].tolist()
  df_result.index = np.random.permutation(range(len(df_result)))
  df_result.sort_index(inplace=True)
  df_result = df_result.sort_values(cn.REPLICATION)
  df_result[cn.GROUP] = groups
  #
  return df_result

def getLines():
  return makeCultureIsolateMutationDF()[
        cn.LINE].unique().tolist()

def stripMutationGroup(mutation_group):
  """
  Mutation groups are specified by a GROUP_STRING.
  Replaces a group with the first mutation in the group.
  :param list-str mutation_group:
  :return list-str:
  """
  results = []
  for element in mutation_group:
    new_element = element.split(MUTATION_GROUP_STRING)[0]
    results.append(new_element)
  return results

def makeIsolateData(**kwargs):
  """
  Creates a dataframe of isolates from genotype_phenotype.
  :return pd.DataFrame:
  """
  df = makeCultureIsolateMutationDF(**kwargs)
  return filterOutlierCultures(df)
