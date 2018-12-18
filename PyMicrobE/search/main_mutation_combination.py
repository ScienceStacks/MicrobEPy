"""Creates files for mutation combinations."""

import __init__
import constants as cn
import pandas as pd
from mutation_combination import MutationCombination
from mutation_context import nextMutationContext

import os

WRITE_COUNT = 5  # Number of combinations between writes
COMBINATION_SIZE_ONE_LINE = 2
COMBINATION_SIZE_ALL_LINE = 2
MIN_IDX = 0
MAX_IDX = 1
IS_MAIN = __name__ == '__main__'


############################################################
# Obtain combinations
############################################################
def doCombination(max_combination, mutation_context,
    lines, out_files, combination_class=MutationCombination):
  """
  Processes a maximum number of combinations, incorporating
  previous results if they are present.
  :param int max_combination: Maximum number of combinations to process
  :param MutationContext mutation_context:
  :param list-str lines:
  :param list-str out_files: output file containing previous results (if they exist)
  :param Type combination_class: class used to process the request
  :return pd.DataFrame, pd.DataFrame: df_min, df_max
    returns an empty dataframe if no new results are present
  """
  def getExistingData(out_file):
    """
    Obtains existing results, if any.
    :param str out_file: file in which data are contained
    :return pd.DataFrame:
    """
    if os.path.isfile(out_file):
      df = pd.read_csv(out_file)
      sel = [r[cn.LINE] in lines for _, r in df.iterrows()]
      df_result = df[sel]
    else:
      df_result = pd.DataFrame()
    return df_result
  #
  # Initialize the min and max dataframes
  df_min_full = getExistingData(out_files[MIN_IDX])
  df_max_full = getExistingData(out_files[MAX_IDX])
  # Find the combinations completed
  if len(df_min_full) > 0:
    combinations = [eval(m) for m in 
        df_min_full[cn.MUTATIONS].unique()]
  else:
    combinations = []
  # Process the request
  combination = combination_class(mutation_context)
  df_min, df_max = combination.do(max_combination, is_tstat=True, 
    is_resample=True, lines=lines, excludes=combinations,
    num_combinations=WRITE_COUNT)
  if len(df_min) != 0:
    df_min_full = pd.concat([df_min_full, df_min])
    df_max_full = pd.concat([df_max_full, df_max])
  else:
    # Signal that everything has been processed
    df_min_full = df_min
    df_max_full = df_max
  return df_min_full, df_max_full


def main(combination_class=MutationCombination):
  """
  Driver for processing mutation combinations.
  Writes results to CSV files.
  :param Type combination_class: Class that processes requests.
  """
  prefix = combination_class.getFilePrefix()
  for ctx in nextMutationContext():
    # Iterate for all combinations of dependent variables and mutation columns
    min_out_file = "%smutation_combination_%s_min.csv" % (prefix, str(ctx))
    max_out_file = "%smutation_combination_%s_max.csv" % (prefix, str(ctx))
    out_files = [min_out_file, max_out_file]  # MIN_IDX, MAX_IDX
    if IS_MAIN:
      print ("Processing %s, %s ..." % (ctx.depvar, ctx.mutation_column))
    #
    def doIterationsForLines(max_size, lines):
      """
      Accumulate results for a single line collection.
      """
      done = False
      while not done:
        df_min, df_max = doCombination(max_size, ctx,
            lines, out_files,
            combination_class=combination_class)
        if len(df_min) == 0:
          done = True
        else:
          df_min.to_csv(out_files[MIN_IDX], index=False)
          df_max.to_csv(out_files[MAX_IDX], index=False)
    #
    doIterationsForLines(COMBINATION_SIZE_ONE_LINE, 
        [cn.LINE_HA2, cn.LINE_HR2, cn.LINE_UE3])
    doIterationsForLines(COMBINATION_SIZE_ALL_LINE, [cn.LINE_ALL])
    

if IS_MAIN:
  main()
