"""Create and manage groups formed from correlated items."""

import __init__
import constants as cn
from correlation_statistic import SetSignificanceLevel
import correlation_statistic as cs
from genome_correlation import GenomeCorrelation
from partition import Partition
import util

import copy
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from sklearn.cluster import KMeans
import warnings


# Clustering algorithms
KMEANS = 'kMeans'
CATCO = 'cdcca'
MAX_SIGLVL = 'max_siglvl'

CSV_CORRELATION_GROUP = "correlation_group.csv"
OUTPUT_DIRECTORY = util.getDataModelPath(None)
FRAC_DEFAULT = 0.5
# Parameters used to generate the groups by line and species
FRAC_DICT = {
    (cn.LINE_HA2, cn.SPECIES_MIX_DVH): 0.1,
    (cn.LINE_HR2, cn.SPECIES_MIX_DVH): 0.15,
    (cn.LINE_UE3, cn.SPECIES_MIX_DVH): 0.1,
    (cn.LINE_HA2, cn.SPECIES_MIX_MMP): 0.1,
    (cn.LINE_HR2, cn.SPECIES_MIX_MMP): 0.1,
    (cn.LINE_UE3, cn.SPECIES_MIX_MMP): 0.1,
    }
CATCO_DICT = {MAX_SIGLVL: 0.001}

# Column names
SIGLVL = "siglvl"


class CorrelationGroup(object):
  """
  Creates groups by clustering items in a correlation matrix.
  Manages these groups.
  """

  def __init__(self, 
      genome_correlation, 
      is_siglvl=True, 
      is_test=False,
      output_directory=OUTPUT_DIRECTORY):
    """
    :param GenomeCorrelation genome_correlation: 
       GenomeCorrelation object
    :param bool is_siglvl: Interpret values as a significance level;
                           otherwise, as correlations
    :param bool is_test: test invocation
    """
    self._genome_correlation = genome_correlation
    self._df_binary = genome_correlation.df_binary
    self._df_corr = genome_correlation.makeCorrelationDF()
    self.categoricals = self._df_corr.columns.tolist()
    self._is_siglvl = is_siglvl
    self.is_test = is_test
    self._output_directory = output_directory
    self._validate()

  def _validate(self):
    if self._df_binary is not None:
      residual = set(
          self._df_corr.columns).symmetric_difference(self._df_binary)
      if len(residual) > 0:
        raise ValueError("df_corr and df_binary must have the same columns")

  def _makeMatrix(self):
    """
    :return np.array: 2-d of correlation values by row
    """
    SMALL_VALUE = 1e-7
    df = self._df_corr.applymap(lambda x: SMALL_VALUE
        if (util.isNull(x) or np.isclose(x,0) or x < 0) else x)
    if self._is_siglvl:
      df = df.applymap(lambda x: -np.log(max(x, SMALL_VALUE)))
    else:
      df = df.applymap(lambda x: -np.log(x))
    return df.values

  def makeKmeansGroups(self, num_groups):
    """
    Clusters the categoricals so that items that are highly
    correlated are grouped together.
    :param int num_groups:
    :return pd.DataFrame: cn.CATEGORICAL, cn.GROUP
    """
    #
    matrix = self._makeMatrix()  # 2-d array with normalized values
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      kmeans = KMeans(
          n_clusters=num_groups, random_state=0).fit(matrix)
    df = pd.DataFrame({
        cn.CATEGORICAL: self.categoricals,
        cn.GROUP: kmeans.labels_
        })
    return df

  def _makeCatCoGroups(self, cluster_parms):
    """
    Creates groups using the categorical data co-occurance clustering algorithm (CATCO).
    The algorithm finds a local minimia for the objective function, the product
    of the significance levels of the sets in the partition.
    :param pd.DataFrame df_matrix: columns are GGENE_ID, rows are KEY_ISOLATE
    :param dict cluster_parms: 
        MAX_SIGLVL: float - maximum significance level for a group
    :return pd.DataFrame:
        cn.GROUP - int group identifier
        cn.CATEGORICAL
        cn.SIGLVL
    """
    max_sl = cluster_parms[MAX_SIGLVL]
    ssl = SetSignificanceLevel(self._df_binary)  # Computes significance levels
    # Consider several values for initializing the partition based on group size.
    fracs = [0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 0.9]
    best_partition = None
    for frac in fracs:
      # Initialize the partition sets to k-means
      num_items = len(self.categoricals)
      num_groups = int(frac*num_items)
      df_groups = self.makeKmeansGroups(num_groups)
      #
      groups = df_groups.groupby(cn.GROUP).groups
      initial_partition = []
      for group in list(groups.keys()):
        indicies = groups[group].values
        this_set = set(df_groups.loc[indicies, cn.CATEGORICAL].tolist())
        initial_partition.append(this_set)
      partition = Partition(set(df_groups[cn.CATEGORICAL]), partition=initial_partition)
      # Iteratively adjust the partition based on significance level
      while True:
        last_sl = ssl.calcLogSets(partition.sets)
        for cat in partition.elements:
          best_set = cn.NULL_SET  # placing cat as a singleton
          best_sl = ssl.calcSet(best_set.union([cat]))
          for a_set in partition.sets:
            mutation_set = partition.move(cat, a_set, is_update=False)
            if a_set != mutation_set.cur_dst:
              new_dst_sl = ssl.calcSet(mutation_set.new_dst)
              if new_dst_sl < best_sl:
                best_set = mutation_set.cur_dst
                best_sl = new_dst_sl
          # See if there is a reduction in the objective function
          lhs_sl = ssl.calcSet(mutation_set.new_src)*best_sl
          rhs_sl = ssl.calcSet(mutation_set.cur_src)*ssl.calcSet(best_set)
          if lhs_sl < rhs_sl:
            mutation_set = partition.move(cat, best_set, is_update=True)
        new_sl = ssl.calcLogSets(partition.sets)
        if np.isclose(last_sl, new_sl):
          break
      if best_partition is None:
        best_partition = partition
      elif ssl.calcLogSets(partition.sets) < ssl.calcLogSets(best_partition.sets):
        best_partition = partition
    # Create the output dataframe
    ssls = [ssl.calcSet(s) for s in partition.sets]
    pairs = zip(partition.sets, ssls)
    categoricals = []
    groups = []
    values_sl = []
    group = 0
    for a_set, value_sl in pairs:
      if value_sl <= max_sl:
        # Valid group
        groups.extend(np.repeat(group, len(a_set)))
      else:
        groups.extend(range(group, group + len(a_set)))
      categoricals.extend(list(a_set))
      values_sl.extend(np.repeat(value_sl, len(a_set)))
      group = max(groups) + 1
    df_result = pd.DataFrame({
        cn.CATEGORICAL: categoricals,
        cn.GROUP: groups,
        cn.SIGLVL: values_sl})
    return df_result

  @classmethod
  def makeCorrelationGroupCSV(cls, 
      cluster_alg=CATCO,
      cluster_parms= CATCO_DICT,
      filename=CSV_CORRELATION_GROUP, 
      output_directory=OUTPUT_DIRECTORY, min_group_size=2,
      min_isolate=2):
    """
    Creates a CSV file with groups for all lines.
    :param str cluster_alg: 'CATCO' or 'kMeans'
    :param dict cluster_parms: parameter used for clustering.
      KMEANS:
        key is (species, line); value is float,
        the fraction of the number of mutations to use as groups
        in clusering (used if 'kMeans')
      CATCO:
        MAX_SIGLVL: float
    :param str filename: name of the file to output.
        No file is written if filename is None
    :param str output_directory: path where output is written
    :param int min_group_size: Minimum size for a group
    :param int min_isolate: Minimum number of isolates
        in which the GGENE_ID must be present to be groupable
    :return pd.DataFrame: 
        cn.SPECIES, cn.LINE, cn.GGENE_ID, 
        cn.GROUP - int group identifier
        cn.COUNT - count of items in the group
        cn.GENE_DESC
    """
    # Find the lines for which the analysis is computed
    query = """
        select distinct line from genotype 
        where transfer = 152 and species is not Null
    """
    df_query = util.readSQL(query)
    dfs = []
    for line in df_query[cn.LINE]:
      for species in [cn.SPECIES_MIX_DVH, cn.SPECIES_MIX_MMP]:
        constraints = [
            lambda r: r[cn.SPECIES] == species,
            lambda r: r[cn.LINE] == line,
            ]
        gc = GenomeCorrelation(is_siglvl=True, 
            constraints=constraints)
        corr_grp = CorrelationGroup(gc)
        if cluster_alg == KMEANS:
          if not (species, line) in list(cluster_parms.keys()):
            frac = FRAC_DEFAULT
          num_groups = max(int(frac*len(corr_grp.categoricals)), 1)
          df_group = corr_grp.makeKmeansGroups(num_groups)
        elif cluster_alg == CATCO:
          df_group = corr_grp._makeCatCoGroups(cluster_parms)
        else:
          import pdb; pdb.set_trace()
          raise ValueError("Invalid cluster algorithm")
        df_group[cn.LINE] = line
        df_group[cn.SPECIES] = species
        # Eliminate from consideration GGENE_ID that occur only once
        query = '''
          select distinct ggene_id, count(distinct key_isolate) as total
                from genotype 
                where ggene_id is not null and species is not null and is_an_mutation = 0
                  and line = '%s' and species = '%s'
                  and is_low_coverage_isolate = 0
          group by ggene_id
        ''' % (line, species)
        df_cnt = util.readSQL(query)
        groupables = df_cnt[df_cnt[cn.TOTAL] >= min_isolate][cn.GGENE_ID].tolist()
        sel = [g in groupables for g in df_group[cn.CATEGORICAL]]
        df_group = df_group.loc[sel,:]
        dfs.append(df_group)
    df_result = pd.concat(dfs, axis=0)
    df_result.rename(columns={cn.CATEGORICAL: cn.GGENE_ID}, inplace=True)
    # Add the counts
    df_result_sub = df_result[[cn.SPECIES, cn.LINE, cn.GROUP, cn.GGENE_ID]]
    df_count = df_result_sub.groupby([cn.SPECIES, cn.LINE, cn.GROUP]).count()
    df_count.reset_index(inplace=True)
    df_count.rename(columns={cn.GGENE_ID: cn.COUNT}, inplace=True)
    df_result = df_result.merge(
        df_count, on=[cn.SPECIES, cn.LINE, cn.GROUP], how='inner')
    df_result = df_result[df_result[cn.COUNT] >= min_group_size]
    # Add the gene description
    ggene_ids = df_result[cn.GGENE_ID].unique().tolist()
    df_desc = util.makeGeneDescriptionDF(ggene_ids)
    del df_desc[cn.GENE_ID]
    df_result = df_result.merge(df_desc, on=cn.GGENE_ID, how='left')
    # Output
    if filename is not None:
      path = os.path.join(output_directory, filename)
      df_result.to_csv(path, index=False)
    return df_result

  def plotGroups(self):
    """
    Plots SSQs versus group size.
    """
    matrix = self._makeMatrix()
    def calc(num):
      with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        kmeans = KMeans(n_clusters=num, random_state=0).fit(matrix)
      return kmeans.inertia_
    #
    sizes = range(1, len(self.categoricals))
    ssqs = [calc(num) for num in sizes]
    # put the major ticks at the middle of each cell
    plt.bar(sizes, ssqs)
    plt.xlabel("No. Groups")
    plt.ylabel("SSQ")
    if not self.is_test:
      plt.show()

  def makeGroupedDF(self,
      cluster_alg=CATCO,
      cluster_parms= CATCO_DICT):
    """
    Changes the columns and rows of the instantiation DF
    so that correlated columns are adjacent.
    :param str cluster_alg:
    :param object cluster_parms: Specific to the algorithm
    :return pd.dataframe: same columns as the base correlation
    """
    if cluster_alg == KMEANS:
      if cluster_parms is None:
        frac = FRAC_DEFAULT
      else:
        frac = cluster_parms
      num_groups = max(int(frac*len(self.categoricals)), 1)
      df_group = self.makeKmeansGroups(num_groups)
    elif cluster_alg == CATCO:
      df_group = self._makeCatCoGroups(cluster_parms)
      del df_group[cn.SIGLVL]
    else:
      raise ValueError("Invalid cluster algorithm")
    df = df_group.sort_values(by=cn.GROUP)
    columns = df[cn.CATEGORICAL]
    df_result = self._df_corr[columns].copy()  # Order the columns
    df_result[cn.CATEGORICAL] = df_result.index.tolist()
    df_result = df_result.merge(df, on=cn.CATEGORICAL, how='inner')
    del df_result[cn.CATEGORICAL]
    df_result = df_result.sort_values(by=cn.GROUP)
    del df_result[cn.GROUP]
    return df_result

  @classmethod
  def makeCorrelationGroupDF(cls, species, line_row, line_col, df=None):
    """
    Computes the statistics about GGENE_IDs common between the cluster
    groups of two lines.
    :param str species:
    :param str line_row: line whose groups are rows in the resulting DF
    :param str line_col: line whose groups are columns in the resulting DF
    :param pd.DataFrame df: cn.GROUP, cn.GGENE_ID, cn.LINE, cn.SPECIES
        if None, uses CSV_CORRELATION_GROUP
    :return dict-of-pd.DataFrame or None:  returns None if no groups
        dict keys are: 
          cn.COUNT - count of GGENE_ID in common
          cn.FRACTION - fraction of GGENE_ID in common
          cn.SET- set of GGENE_ID in common
          cn.SIGLVL - significance level of the correlation
        The dataframes have the columns:
          index are groups in line_row,
          groups in line_col
    """
    KEYS = [cn.COUNT, cn.FRACTION, cn.SET, cn.SIGLVL]
    def makeGroupDF(line):
      sel = [(r[cn.SPECIES] in species) and (r[cn.LINE] in line)
             for _,r in df.iterrows()]
      df_group = df.loc[sel, :]
      groups = df_group[cn.GROUP].unique()
      groups.sort()
      return df_group, groups
    def makeGgeneSet(df, group):
      return set(df[df[cn.GROUP] == group][cn.GGENE_ID])
    #
    if df is None:
      path = os.path.join(OUTPUT_DIRECTORY, CSV_CORRELATION_GROUP)
      df = pd.read_csv(path)
    # Construct the group dataframe
    df_row, group_rows  = makeGroupDF(line_row)
    df_col, group_cols  = makeGroupDF(line_col)
    # Check for null values
    if (len(group_rows) == 0) or (len(group_cols) == 0):
      return None
    # Compute the group overlaps
    rows_dict = {k: [] for k in KEYS}
    for group_row in group_rows:
      set_row = makeGgeneSet(df_row, group_row)
      row_dict = {k: {cn.INDEX: group_row} for k in KEYS}
      for group_col in group_cols:
        set_col = makeGgeneSet(df_col, group_col)
        denom = min(len(set_row), len(set_col))
        common = set_row.intersection(set_col)
        row_dict[cn.SET][group_col] = common
        row_dict[cn.COUNT][group_col] = len(common)
        row_dict[cn.FRACTION][group_col] = 1.0*len(common) / denom
        row_dict[cn.SIGLVL][group_col] = cls.calcIntergroupSiglvl(
            species, line_row, line_col, group_row, group_col)
      [rows_dict[k].append(row_dict[k]) for k in KEYS]
    dfs = {k: pd.DataFrame(rows_dict[k]) for k in KEYS}
    [dfs[k].set_index('index', inplace=True) for k in KEYS]
    return dfs

  @classmethod
  def calcIntergroupSiglvl(cls, species, line1, line2,
      group1, group2, df=None):
    """
    Calculates the significance level for the overlap of
    GGENE_ID of two groups in different lines.
    :param str species:
    :param str line1:
    :param int group1:
    :param str line2:
    :param int group2:
    :param pd.DataFrame df: cn.GROUP, cn.GGENE_ID, cn.LINE, cn.SPECIES
        if None, uses CSV_CORRELATION_GROUP
    :return float, float: lower and upper bounds on significance level
    """
    def getMutations(df, line, group=None):
      """
      :param pd.DataFrame df: cn.LINE, cn.GROUP, cn.GGENE_ID (opt)
      :param str line:
      :param int group:
      :return set: set of cn.GGENE_ID
      """
      if group is not None:
        sel = [(r[cn.LINE] == line) and (r[cn.GROUP]  == group) 
               for _,r in df.iterrows()]
      else:
        sel = [(r[cn.LINE] == line) for _,r in df.iterrows()]
      result = set(df.loc[sel, cn.GGENE_ID])
      return result
    #
    if df is None:
      path = os.path.join(OUTPUT_DIRECTORY, CSV_CORRELATION_GROUP)
      df = pd.read_csv(path)
    #
    df_base = GenomeCorrelation.makeBaseDF()
    # Limit to the species of interest
    df_base = df_base[df_base[cn.SPECIES] == species]
    # Compute the counts of items in common
    ggenes1 = getMutations(df_base, line1)
    ggenes2 = getMutations(df_base, line2)
    M = len(ggenes1.intersection(ggenes2))
    # Compute counts for the size of each group
    ggenes1 = getMutations(df, line1, group1)
    n1 = len(ggenes1)
    ggenes2 = getMutations(df, line2, group2)
    n2 = len(ggenes2)
    # Compute the size of the intersection of the two groups
    k = len(ggenes1.intersection(ggenes2))
    #
    probs = cs.calcCumlProbIntersection(M, n1, n2, k)
    return probs
        

if __name__ == '__main__':
  print("***Creating the correlation groups.")
  CorrelationGroup.makeCorrelationGroupCSV()
