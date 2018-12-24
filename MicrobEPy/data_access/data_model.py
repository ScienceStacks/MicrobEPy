"""
Creates the Data Model consisting of the following tables
that are written as CSV files.
  culture, isolate, mutation, gene_description
  culture_isolate_link, isolate_mutation_link,
Details of the columns in the tables are found in the constants
file in the definitions of TABLE_SCHEMAS.

The key for culture, isolate, mutation are, respectively,
  KEY_CULTURE  - Constructed. "C" + an integer
  KEY_ISOLATE  - Composed of multiple elements, as
                 described in isolate.py
  KEY_MUTATION - gene_id.position.initial_nt.changed_nt
    initial_nt is required because mutations at the same position
    can be of different lengths.
  GENE_ID - a key for gene_description

Navigation between the files is done by the above keys. For example,
we navigate from culture to culture_isolate_link via KEY_CULTURE.
The "link" files provide a way to join culture with isolate, mutation
with isolate, and so join mutation with culture.

Columns in these files are obtained either from growth or sequence
data, or are constructured in this module.

Data sets produced
  mutation.csv - all columns pertaining to a mutation
  isolate.csv - all columns relating to a community (e.g., isolate, EPD, line)
  culture.csv - all columns relating to a culture that is assayed for growth and yield
  isolate_mutation_link.csv - relates isolates to mutations and includes information
     on mutation fractions for isolates
  culture_isolate_link.csv - relates cultures to isolates
  culture_isolate_mutation.csv - "outer" join of culture, isolate, mutation. Has "null"
     values where a value is missing (e.g., the rate values for an isolate that was not cultured)
  gene_description.csv - Information on genes
  genotype.csv - mutation + isolate + gene information
  genotype_phenotype.csv - "inner" join of mutation, isolate, and culture (via the link columns)
     excludes isolates and cultures for which there is not sequence data (e.g., low coverage isolates)

All of the foregoing also appear as tables in convolution.db

"""

import __init__
import constants as cn
from isolate import Isolate
import isolate as iso
import util
import util_data_access

import os
import pandas as pd
import numpy as np


DVH_MUTATION_FILE = "dvh_mutations_allsamples_attributes_09042018.xlsx"
MMP_MUTATION_FILE = "mmp_mutations_allsamples_attributes_09042018.xlsx"
DVH_PAIRING_FILE = "names_clonal_isolates_dvh.csv"
MMP_PAIRING_FILE = "names_clonal_isolates_mmp.csv"
# Columns: rate, yield, isolate_1, and optionally isolate_2
RATE_YIELD_FILES = ["early_evolved_lines.xlsx",
    "epds_growth.xlsx",
    "growth_data.xlsx",
    "growth_data_wt_pairings.xlsx",
     ]

# Columns
KEY_ISOLATE_FRONT = "key_isolate_front"
KEY_ISOLATE_BACK = "key_isolate_back"
FREQ1 = "freq1"
FREQ2 = "freq2"
OTHER_PREDICTOR = "other_predictor"
MUTATION_DELETES = [FREQ1, FREQ2, OTHER_PREDICTOR]

# Tables
TABLE_ISOLATE_MUTATION = "isolate_mutation"
TABLE_CULTURE_ISOLATE = "culture_isolate"

###########################################################
# HELPER FUNCTIONS
############################################################
def makeKeyMutation(df):
  """
  The mutation key consists of the components
    GENE_ID, POSITION, INITIAL_NT, CHANGED_NT
  Creates values for KEY_MUTATION, which is
  gene (or DVH/MMP for intergenic)_POSITION_CHANGED_NT
  :param pd.DataFrame df: cn.GENE_ID, cn.POSITION, cn.CHANGED_NT
  :return list-of-str:
  """
  tuples = zip(df[cn.GENE_ID].tolist(),
      df[cn.POSITION].tolist(),
      df[cn.INITIAL_NT].tolist(),
      df[cn.CHANGED_NT].tolist(),
      df[cn.SPECIES].tolist())
  format_dict = {cn.SPECIES_MIX_DVH: 'DVH__',
       cn.SPECIES_MIX_MMP: 'MMP__'}
  result = []
  for gene_id, position, initial_nt, changed_nt, species in tuples:
    if util.isNull(gene_id):
      gene_id = format_dict[species]
    result.append("%s%s%s%s%s%s%s" % (
        gene_id, cn.MUTATION_SEPARATOR, str(position),
        cn.MUTATION_SEPARATOR, initial_nt,
        cn.MUTATION_SEPARATOR, changed_nt))
  return result

def makeRawPairingDF(csv_file):
  """
  :param str filename: pairing file
  :return pd.DataFrame: cn.SAMPLE, KEY_ISOLATE_FRONT, KEY_ISOLATE_BACK
  """
  df = pd.read_csv(util.getSequenceDataPath(csv_file),
     sep=",")
  current_columns = df.columns.tolist()
  new_columns = [cn.SAMPLE, KEY_ISOLATE_FRONT, KEY_ISOLATE_BACK]
  for nn in range(len(current_columns)):
    df.rename(columns={current_columns[nn]: new_columns[nn]},
        inplace=True)
  return df

def makeSpeciesMixInCulture(df):
  """
  :param pd.DataFrame df: has KEY_ISOLATE, KEY_CULTURE
  :return pd.DataFrame: Columns are: cn.KEY_CULTURE, cn.SPECIES_MIX
  """
  # Dictionary maps the isolates present to the species mix
  MONOCULTURE = "monoculture"
  COMMUNITY = {
      (cn.ISOLATE_CLONE, cn.ISOLATE_CLONE): cn.SPECIES_MIX_BOTH,
      (cn.ISOLATE_EPD, cn.ISOLATE_UNKNOWN): cn.SPECIES_MIX_EPD,
      (cn.ISOLATE_LINE, cn.ISOLATE_UNKNOWN): cn.SPECIES_MIX_LINE,
      (cn.ISOLATE_LINE, cn.ISOLATE_LINE): cn.SPECIES_MIX_LINE,
      (cn.ISOLATE_CLONE, cn.ISOLATE_UNKNOWN): MONOCULTURE,
      (cn.ISOLATE_CLONE, cn.ISOLATE_LINE): cn.SPECIES_MIX_BOTH,  # Wild type pairing
      (cn.ISOLATE_UNKNOWN, cn.ISOLATE_UNKNOWN): cn.SPECIES_MIX_UNKNOWN,
      }
  # Construct the mixes
  mixes = []
  for stg1, stg2 in zip(df[cn.ISOLATE1], df[cn.ISOLATE2]):
    isolate1 = Isolate.create(stg1)
    isolate2 = Isolate.create(stg2)
    pair = (isolate1.getCommunity(), isolate2.getCommunity())
    mix = COMMUNITY[pair]
    if mix == MONOCULTURE:
      mix = isolate1.species
    mixes.append(mix)
  return pd.DataFrame({
      cn.KEY_CULTURE: df[cn.KEY_CULTURE],
      cn.SPECIES_MIX: mixes
      })

def deleteUnneededMutationColumns(df):
  """
  Deletes columns from the mutation dataframe that are not used.
  :param pd.DataFrame df: mutation dataframe
  """
  for col in MUTATION_DELETES:
    if col in df.columns:
      del df[col]


############################################################
# Data Model
############################################################
class DataModel(object):

  def __init__(self,
      dvh_mutation_file=DVH_MUTATION_FILE,
      mmp_mutation_file=MMP_MUTATION_FILE,
      dvh_pairing_file=DVH_PAIRING_FILE,
      mmp_pairing_file=MMP_PAIRING_FILE,
      rate_yield_files=RATE_YIELD_FILES,
      output_directory=None,
      ):
    if output_directory is None:
      output_directory = util.getDataModelPath(None)
    self._dvh_mutation_file = dvh_mutation_file
    self._mmp_mutation_file = mmp_mutation_file
    self._dvh_pairing_file = dvh_pairing_file
    self._mmp_pairing_file = mmp_pairing_file
    self._rate_yield_files = rate_yield_files
    self._output_directory = output_directory
    self._dfs = None

  def do(self, is_print=False):
    """
    Create the CSV files for the data model.
    :param bool is_print: flag to print statistics
    """
    self._makeCombinationDFs()
    # Create the data model
    self._makeBaseDataModelDFs()
    self._makeJoinedDataModelDFs()
    # Save the results
    self._writeDFs(is_print=is_print)

  def _makeCombinationDFs(self):
    """
    These data frames are modest transformations of the input data.
    """
    self._dfs = {}
    self._dfs[TABLE_CULTURE_ISOLATE] = self._makeCultureIsolateDF()
    self._dfs[TABLE_ISOLATE_MUTATION] = self._makeIsolateMutationDF()

  def _makeGenotypePhenotypeDF(self):
    """
    Constructs the merge of mutation, isolate, culture, and gene_description.
    Only those isolates that have mutation information AND are cultured
    are included.
    Only cultures that have sequenced isolates are included.
    """
    columns = cn.TABLE_SCHEMAS.getColumns([cn.TABLE_GENOTYPE_PHENOTYPE])
    df = self._dfs[cn.TABLE_CULTURE_ISOLATE_MUTATION][columns].copy()
    constraint = lambda r: (not util.isNull(r[cn.KEY_MUTATION]))  \
        and (not util.isNull(r[cn.KEY_ISOLATE]))  \
        and (not util.isNull(r[cn.KEY_CULTURE]))  \
        and (not r[cn.KEY_CULTURE] in cn.LOW_COVERAGE_ISOLATES)
    sel = df.apply(constraint, axis=1)
    df_result = df.loc[sel].copy(deep=True)
    df_result.drop_duplicates(inplace=True)
    return df_result

  def _makeGenotypeDF(self):
    """
    Constructs the merge of mutation, isolate, and gene_description.
    """
    columns = cn.TABLE_SCHEMAS.getColumns(
        [cn.TABLE_GENOTYPE])
    df = self._dfs[cn.TABLE_CULTURE_ISOLATE_MUTATION][columns]
    constraint = (lambda r: 
        (not util.isNull(r[cn.KEY_ISOLATE])) and
        (not util.isNull(r[cn.KEY_MUTATION]))
        )
    sel = df.apply(constraint, axis=1)
    df_result = df.loc[sel].copy(deep=True)
    df_result.drop_duplicates(inplace=True)
    return df_result

  def _makeBaseDataModelDFs(self):
    """
    Creates data frames dervied from the combination data frames.
    These methods only depend on the combination dataframes
    TABLE_CULTURE_ISOLATE, TABLE_ISOLATE_MUTATION
    """
    self._dfs[cn.TABLE_MUTATION] = self._makeMutationDF()
    self._dfs[cn.TABLE_CULTURE] = self._makeCultureDF()
    self._dfs[cn.TABLE_ISOLATE] = self._makeIsolateDF()
    self._dfs[cn.TABLE_CULTURE_ISOLATE_LINK] = self._makeCultureIsolateLinkDF()
    self._dfs[cn.TABLE_ISOLATE_MUTATION_LINK] = self._makeIsolateMutationLinkDF()
    self._dfs[cn.TABLE_GENE_DESCRIPTION] = self._makeGeneDescriptionDF()

  def _makeJoinedDataModelDFs(self):
    self._dfs[cn.TABLE_CULTURE_ISOLATE_MUTATION] =  \
        self._makeCultureIsolateMutationDF()
    self._dfs[cn.TABLE_GENOTYPE] = self._makeGenotypeDF()
    self._dfs[cn.TABLE_GENOTYPE_PHENOTYPE] =  \
        self._makeGenotypePhenotypeDF()

  def _writeDFs(self, is_print=False):
    """
    Writes the list of data frames as CSVs and as tables in a database.
    """
    #
    if is_print:
      print("DF\t\t\t\tsize(KB)")
    total = 0
    conn = util.getDBConnection()
    for table_name in list(self._dfs.keys()):
      df = self._dfs[table_name]
      kilobytes = int(self._dfs[table_name].memory_usage().sum()/1e3)
      total += kilobytes
      tabs = "\t"
      if len(table_name) < 22:
        tabs += "\t"
      if len(table_name) < 14:
        tabs += "\t"
      if is_print:
        print("  %s%s%d" % (table_name, tabs, kilobytes))
      self._write(df, util_data_access.makeCSVFilename(table_name))
      try:
        df.to_sql(table_name, conn, if_exists="replace",
            index=False)
      except Exception as err:
        import pdb; pdb.set_trace()
        print (err)
    conn.close()
    if is_print:
      print("Total storage(KB):\t\t%d" % total)

  def _makeCultureIsolateMutationDF(self):
    """
    Constructs the join of the culture, isolate, and mutation data
    along with GENE_DESCRIPTION.
    """
    def trimIsolateData(df):
      isolates = [Isolate.create(s) for s in df[cn.KEY_ISOLATE]]
      df[cn.LINE] = [cn.NONE if i.line == cn.ISOLATE_DEFAULT
          else i.line for i in isolates]
      df[cn.TRANSFER] = [cn.NONE if i.transfer == cn.ISOLATE_DEFAULT
          else i.transfer for i in isolates]
      df[cn.EPD_ID] = [cn.NONE if i.epd_id == cn.ISOLATE_DEFAULT
          else i.epd_id for i in isolates]
      df[cn.CLONE] = [cn.NONE if i.clone == cn.ISOLATE_DEFAULT
          else i.clone for i in isolates]
      df[cn.SPECIES] = [cn.NONE if i.species == cn.ISOLATE_DEFAULT
          else i.species for i in isolates]
      df[cn.EXPERIMENT] = [cn.NONE if i.experiment == cn.ISOLATE_DEFAULT
          else i.experiment for i in isolates]
    # Because the universal table is so large, it is assembled in
    # 3 pieces. Two peices are done in-memory. The third piece is a
    # View.
    # Merge one set
    df_result1 = self._dfs[cn.TABLE_MUTATION].merge(
        self._dfs[cn.TABLE_ISOLATE_MUTATION_LINK],
        on=cn.KEY_MUTATION, how='outer')
    # Merge another set
    df_result2 = self._dfs[cn.TABLE_ISOLATE].copy()
    trimIsolateData(df_result2)
    df_result2 = df_result2.merge(
        self._dfs[cn.TABLE_CULTURE_ISOLATE_LINK],
        on=cn.KEY_ISOLATE, how='outer')
    df_result2 = df_result2.merge(self._dfs[cn.TABLE_CULTURE],
        on=cn.KEY_CULTURE, how='outer')
    # Merge the two parts and clean them up.
    df_result = df_result1.merge(df_result2, on=cn.KEY_ISOLATE, how='outer')
    df_result = util.pruneNullRows(df_result)
    df_result = util.mergeLeft(df_result, self._dfs[cn.TABLE_GENE_DESCRIPTION], cn.GENE_ID)
    df_result.drop_duplicates(inplace=True)
    return df_result

  def _renameGene(self, df):
    df[cn.GENE_ID] = df[cn.GENE_ID].apply(
      lambda v: v if util.isNull(v) else v.replace("_", ""))
    return df

  def _makeGeneDescriptionDF(self):
    """
    Creates a dataframe with gene descriptions.
    :return pd.DataFrame:
        see cn.TABLE_SCHEMAS for cn.TABLE_GENE_DESCRIPTION
    """
    def getGeneDataFrame(filename):
      df = util.makeDataframeFromXlsx(
          util.getReferenceDataPath(filename))
      df.rename(columns={
          'gene_ID': cn.GENE_ID,
          'locusId': cn.LOCUS_ID,
          'accession': cn.ACCESSION,
          'GI': cn.GENE_GI,
          'scaffoldId': cn.SCAFFOLD_ID,
          'start': cn.START,
          'stop': cn.STOP,
          'strand': cn.STRAND,
          'name': cn.GENE_NAME,
          'desc': cn.GENE_DESC,
          'COG': cn.COG,
          'COGDesc': cn.COG_DESC,
          'COGFun': cn.COG_FUN,
          'TIGRFam': cn.TIGR_FAM,
          'TIGRRoles': cn.TIGR_ROLES,
          'GO': cn.GO,
          'EC': cn.EC,
          'ECDesc': cn.EC_DESC,
          }, inplace=True)
      self._renameGene(df)
      return df
    #
    df_dvh = getGeneDataFrame("dvh_gene_description.xlsx")
    df_mmp = getGeneDataFrame("mmp_gene_description.xlsx")
    df_result = pd.DataFrame(df_dvh)
    df_result = df_result.append(df_mmp)
    util.replaceNan(df_result, value=None)
    util.cleanDF(df_result, is_reset_index=True)
    return df_result

  def _makeMutationDF(self):
    """
    Constructs the Mutation dataframe by trimming the
    isolate_mutation dataframe and adding a column for 
    indicating ancestral mutations.
    :return pd.DataFrame: df_mutation
    """
    df = self._dfs[TABLE_ISOLATE_MUTATION]
    # Handle AN mutation. A mutation is ancestral if
    # there is a wild type with that mutation.
    an_mutations = [r[cn.KEY_MUTATION]
        for _, r in df.iterrows()
        if Isolate.isAN(r[cn.KEY_ISOLATE])]
    df[cn.IS_AN_MUTATION] = [r[cn.KEY_MUTATION] in an_mutations
        for _, r in df.iterrows()]
    #
    keep_columns = list(cn.TABLE_SCHEMAS.getSchema(
        cn.TABLE_MUTATION).columns)
    df_result = util.trimDF(df, keep_columns=keep_columns)
    #
    util.cleanDF(df_result, is_reset_index=True)
    return df_result

  def _makeCultureDF(self):
    """
    :return pd.DataFrame: df_culture
    """
    df = self._dfs[TABLE_CULTURE_ISOLATE]
    df_result = util.trimDF(df,
        keep_columns=[cn.RATE, cn.YIELD, cn.KEY_CULTURE,
        cn.ISOLATE1, cn.ISOLATE2])
    df_mix = makeSpeciesMixInCulture(df_result)
    df_result = df_result.merge(
        df_mix, on=cn.KEY_CULTURE, how='inner')
    del df_result[cn.ISOLATE1]
    del df_result[cn.ISOLATE2]
    util.cleanDF(df_result, is_reset_index=True)
    return df_result

  def _write(self, df, csv_file):
    """
    Writes a dataframe as a csv file.
    :param pd.DataFrame df:
    :param str csv_file:
    """
    path = os.path.join(self._output_directory, csv_file)
    df.to_csv(path, index=False)

  def _makeIsolateMutationDF(self):
    """
    Creates a single dataframe for the mutation files.
    Adds the columns: cn.KEY_MUTATION, cn.SPECIES, cn.KEY_ISOLATE
    :return pd.DataFrame:
    """
    df_dvh = util.makeDataframeFromXlsx(
        util.getSequenceDataPath(self._dvh_mutation_file))
    deleteUnneededMutationColumns(df_dvh)
    self._renameGene(df_dvh)
    df_dvh[cn.SPECIES] = cn.SPECIES_MIX_DVH
    df_mmp = util.makeDataframeFromXlsx(
        util.getSequenceDataPath(self._mmp_mutation_file))
    deleteUnneededMutationColumns(df_mmp)
    self._renameGene(df_mmp)
    df_mmp[cn.SPECIES] = cn.SPECIES_MIX_MMP
    df = pd.concat([df_dvh, df_mmp], sort=True)
    df[cn.KEY_MUTATION] = makeKeyMutation(df)
    df[cn.KEY_ISOLATE] = self._makeIsolatesForMutations(df)
    df[cn.GENE_POSITION] = ['.'.join(k.split('.')[0:2]) 
        for k in df[cn.KEY_MUTATION]]
    df[cn.GGENE_ID] = [
        r[cn.GENE_ID] if not util.isNull(r[cn.GENE_ID])
        else r[cn.GENE_POSITION]
        for _,r in df.iterrows()]
    return df

  def _makeCultureIsolateDF(self):
    """
    Creates a single dataframe for the rate-yield files.
    :return pd.DataFrame: cn.RATE, cn.YIELD, cn.KEY_CULTURE,
        cn.ISOLATE1, cn.ISOLATE2
    """
    global culture_index
    culture_index = 0
    #
    def createIsolate(series):
      return [s if util.isNull(s) else str(Isolate.create(s))
              for s in series]
    #
    def makeCultureKey():
      global culture_index
      culture_index += 1
      return "C%d" % culture_index
    #
    dfs = []
    for gfile in self._rate_yield_files:
      df = util.makeDataframeFromXlsx(util.getRateYieldDataPath(gfile))
      df[cn.ISOLATE1] = createIsolate(df[cn.ISOLATE1])
      if cn.ISOLATE2 in df.columns:
        df[cn.ISOLATE2] = createIsolate(df[cn.ISOLATE2])
      df[cn.KEY_CULTURE] =  \
          [makeCultureKey() for _ in range(len(df.index))]
      dfs.append(df)
    return pd.concat(dfs, sort=True)

  def _getRawIsolatesInRateYieldFiles(self):
    """
    :return list-of-str: list of isolates in the rate-yield files.
    """
    raw_isolates = []
    for gfile in self._rate_yield_files:
      df = util.makeDataframeFromXlsx(util.getRateYieldDataPath(gfile))
      raw_isolates.append(df[cn.ISOLATE1].tolist())
      if ISOLATE2 in df.columns:
        raw_isolates.append(df[cn.ISOLATE2].tolist())
    return raw_isolates

  def _makeIsolatesForMutations(self, df):
    """
    Creates the isolates for a mutation dataframe.
    :param pd.DataFrame df_mutation: cn.SAMPLE, cn.SPECIES
    :return list-of-str: isolate strings
    """
    isolate_stgs = []
    pairs = zip(df[cn.SAMPLE], df[cn.SPECIES])
    df_pairing = self._makePairingDF()
    for sample, species in pairs:
      clone = cn.ISOLATE_DEFAULT
      experiment = cn.ISOLATE_DEFAULT
      split = sample.split("_")
      # AN_Coculture-Ancestor
      if sample == "AN_Coculture-Ancestor":
        line = cn.LINE_ANC
        transfer = cn.ISOLATE_DEFAULT
        epd_id = cn.ISOLATE_DEFAULT
        species = cn.ISOLATE_DEFAULT
      # Ex: AN_Dv-Ancestor-1
      elif "Ancestor-" in sample:
        if "-1" in sample:
          line = cn.LINE_AN1
        elif "-2" in sample:
          line = cn.LINE_AN2
        else:
          raise ValueError("Invalid AN_Ancestor %s" % sample)
        transfer = cn.ISOLATE_DEFAULT
        epd_id = cn.ISOLATE_DEFAULT
        if "Dv" in sample:
          species = cn.SPECIES_MIX_DVH
        elif "Mm" in sample:
          species = cn.SPECIES_MIX_MMP
        else:
          raise ValueError("Did not find species in ancestor sample: %s"
              % sample)
      # Ex: "EP_HA2_01"
      elif (len(split) == 3) and split[0] == "EP":
        line = split[1]
        transfer = cn.TRANSFER_1000G
        epd_id = split[2]
        species = cn.ISOLATE_DEFAULT
      # Ex: "EP_UE3"
      elif (len(split) == 2) and split[0] == "EP":
        line = split[1]
        transfer = cn.TRANSFER_1000G
        epd_id = "00"
        species = cn.ISOLATE_DEFAULT
      # Ex: "TG_UE3"
      elif (len(split) == 2) and split[0] == "TG":
        transfer = cn.TRANSFER_1000G
        line = split[1]
        epd_id = cn.ISOLATE_DEFAULT
        species = cn.ISOLATE_DEFAULT
      # Ex: "EG_HR2-45"
      elif (len(split) == 2) and split[0] == "EG":
        sub_split = split[1].split("-")
        if len(sub_split) != 2:
          raise ValueError("%s is an invalid sample format"
              % sample)
        line = sub_split[0]
        transfer = sub_split[1]
        epd_id = cn.ISOLATE_DEFAULT
        species = cn.ISOLATE_DEFAULT
      # Ex: "EG_HR2_45"
      elif (len(split) == 3) and (split[0] == "EG"):
        line = split[1]
        transfer = split[2]
        epd_id = cn.ISOLATE_DEFAULT
        species = cn.ISOLATE_DEFAULT
      # Ex: "CI_20_S76"
      elif (len(split) == 3) and split[0] == "CI":
        sel = (df_pairing[cn.SAMPLE] == int(split[1])) &  \
            (df_pairing[cn.SPECIES] == species)
        df1 = df_pairing.loc[sel]
        if len(df1.index) == 0:
          msg = "***Cannot match sample number with a pairing."
          msg += "\nSample is %s and species is %s\n"   \
          % (sample, species)
          print(msg)
          line = cn.ISOLATE_DEFAULT
          transfer = cn.ISOLATE_DEFAULT
          epd_id = cn.ISOLATE_DEFAULT
          clone = cn.ISOLATE_DEFAULT
          experiment = cn.EXPERIMENT_CI
        elif len(df1.index) > 1:
          import pdb; pdb.set_trace()
          msg = "%s: should not have than one row for a match"  \
              % sample
          raise RuntimeError(msg)
        else:
          isolate = Isolate.create(df1[cn.KEY_ISOLATE].tolist()[0])
          line = isolate.line
          transfer = isolate.transfer
          epd_id = isolate.epd_id
          clone = isolate.clone
          species = isolate.species
          experiment = isolate.experiment
      # Ex: "EP_00_S76" - ignore these
      elif (len(split) == 3) and (split[0] == "EP")  \
          and (split[1] == "00"):
        import pdb; pdb.set_trace()
        raise RuntimeError("%s: Should not have '00' samples" %
            sample)
      else:
        import pdb; pdb.set_trace()
        raise ValueError("%s is an invalid sample format"
            % sample)
      isolate_stg = str(Isolate(line=line, transfer=transfer,
          epd_id=epd_id, clone=clone, species=species,
          experiment=experiment))
      isolate_stgs.append(isolate_stg)
    return isolate_stgs

  def _makePairingDF(self):
    """
    :return pd.DataFrame: cn.KEY_ISOLATE, cn.SAMPLE, cn.SPECIES
    """
    def process(pairing_file, species, transfer=cn.TRANSFER_DEFAULT):
      """
      :param str pairing_file:
      :param str species:
      :param int transfer:
      :return pd.DataFrame: cn.SAMPLE, cn.KEY_ISOLATE
      """
      isolates = []
      df_pairing = makeRawPairingDF(pairing_file)
      isolate_fronts = df_pairing[KEY_ISOLATE_FRONT].tolist()
      isolate_backs = df_pairing[KEY_ISOLATE_BACK].tolist()
      pairs = zip(isolate_fronts, isolate_backs)
      for front, back in pairs:
        if front == cn.LINE_AN:
          isolate_stg = "%s%s%s" % (cn.LINE_AN, cn.ISOLATE_SEPARATOR,
              back)
        else:
          front_split = front.split(cn.ISOLATE_SEPARATOR)
          if len(front_split) != 2:
            raise ValueError("%s has wrong format in pairing file."
                % front_split)
          line = front_split[0]
          epd_id = front_split[1]
          isolate_stg = "%s%s%s%s%s%s%s"  \
              % (line, cn.ISOLATE_SEPARATOR,
              str(transfer), cn.ISOLATE_SEPARATOR,
              epd_id, cn.ISOLATE_SEPARATOR, back)
        isolates.append(isolate_stg)
      df_pairing[cn.KEY_ISOLATE] =  \
          [str(Isolate.create(s)) for s in isolates]
      df_pairing[cn.SPECIES] = species
      return df_pairing

    df_dvh = process(self._dvh_pairing_file,
         cn.SPECIES_MIX_DVH)
    df_mmp = process(self._mmp_pairing_file,
         cn.SPECIES_MIX_MMP)
    df = pd.concat([df_dvh, df_mmp], sort=True)
    del df[KEY_ISOLATE_FRONT]
    del df[KEY_ISOLATE_BACK]
    return df

  def _makeIsolateDF(self):
    """
    Creates a single dataframe with the isolates in the standard
    isolate format, along with related columns.
    :return pd.DataFrame:
        See TABLE_SCHEMAS for TABLE_ISOLATE
    Notes:
        1. '*' values are converted to None so that they appear as
           empty fields in the CSV file.
    """
    isolate_stgs = []
    for table in [TABLE_ISOLATE_MUTATION, TABLE_CULTURE_ISOLATE]:
      df = self._dfs[table]
      if cn.KEY_ISOLATE in df.columns:
        isolate_stgs.extend(df[cn.KEY_ISOLATE].unique().tolist())
      if cn.ISOLATE1 in df.columns:
        isolate_stgs.extend(df[cn.ISOLATE1].unique().tolist())
      if cn.ISOLATE2 in df.columns:
        isolate_stgs.extend(df[cn.ISOLATE2].unique().tolist())
    isolates = [Isolate.create(s) for s in set(isolate_stgs)
                if not util.isNull(s)]
    df_result = pd.DataFrame({
        cn.KEY_ISOLATE: [str(i) for i in isolates],
        cn.LINE: [i.line for i in isolates],
        cn.TRANSFER: [i.transfer for i in isolates],
        cn.EPD_ID: [i.epd_id for i in isolates],
        cn.CLONE: [i.clone for i in isolates],
        cn.SPECIES: [i.species for i in isolates],
        cn.EXPERIMENT: [i.experiment for i in isolates],
        cn.IS_LOW_COVERAGE_ISOLATE: False,
        })
    for col in [cn.LINE, cn.TRANSFER, cn.EPD_ID,
        cn.CLONE, cn.SPECIES, cn.EXPERIMENT]:
      df_result[col] = df_result[col].apply(
          lambda v: None if v==cn.ISOLATE_DEFAULT else v)
    sel = df_result[cn.KEY_ISOLATE].isin(cn.LOW_COVERAGE_ISOLATES)
    df_result.loc[sel, cn.IS_LOW_COVERAGE_ISOLATE] = True
    util.cleanDF(df_result, is_reset_index=True)
    return df_result

  def _makeCultureIsolateLinkDF(self):
    """
    Create a DataFrame that contains information related to
    the combination of a culture and isolate.
    Uses TABLE_CULTURE_ISOLATE which has the columns
        cn.ISOLATE1, cn.ISOLATE2, cn.KEY_CULTURE
    :return pd.DataFrame:
    """
    df_culture_isolate = self._dfs[TABLE_CULTURE_ISOLATE]
    df1 = util.trimDF(df_culture_isolate,
        [cn.KEY_CULTURE, cn.ISOLATE1])
    df1.rename(columns={cn.ISOLATE1: cn.KEY_ISOLATE}, inplace=True)
    df2 = util.trimDF(df_culture_isolate,
        [cn.KEY_CULTURE, cn.ISOLATE2])
    df2 = pd.DataFrame([r for _,r in df2.iterrows()
                        if not util.isNull(r[cn.ISOLATE2])])
    df2.rename(columns={cn.ISOLATE2: cn.KEY_ISOLATE}, inplace=True)
    return pd.concat([df1, df2], sort=True)

  def _makeIsolateMutationLinkDF(self):
    """
    Create a DataFrame that contains information related to
    the combination of an isolate and a mutation.
    :return pd.DataFrame:
    """
    def floatFormat(stg):
      return float(stg.replace("%",""))
    #
    SEP = ':'
    PREDICTORS = set([cn.GATK, cn.SAMTOOLS, cn.VARSCAN])
    #
    df_culture_isolate = self._dfs[TABLE_CULTURE_ISOLATE]
    df_isolate_mutation = self._dfs[TABLE_ISOLATE_MUTATION]
    #
    predictor_dict = {p: [] for p in PREDICTORS}
    keep_columns = [cn.KEY_ISOLATE, cn.KEY_MUTATION, cn.READ_NUMBER,
        cn.PREDICTOR, cn.FREQ_PREDICTOR, cn.FREQ]
    df = util.trimDF(df_isolate_mutation, keep_columns=keep_columns)
    #
    df[cn.FREQ] = df[cn.FREQ].apply(floatFormat)
    # Construct to predictor columns
    pairs = zip(df[cn.PREDICTOR].tolist(),
        df[cn.FREQ_PREDICTOR].tolist())
    for predictor, freq_predictor in pairs:
      predictor_split = predictor.split(SEP)
      freq_predictor_split = freq_predictor.split(SEP)
      if len(predictor_split) != len(freq_predictor_split):
        msg = "%s, %s: len(predictor) != len(freq_predictor"  \
            % (predictor, freq_predictor)
        raise RuntimeError(msg)
      # Handle predictors present
      for prd, frq in zip(predictor_split, freq_predictor_split):
        val = floatFormat(frq)
        predictor_dict[prd].append(val)
      # Handle missing predictors
      for prd in PREDICTORS.difference(predictor_split):
        predictor_dict[prd].append(0)
    for predictor in PREDICTORS:
      df[predictor] = predictor_dict[predictor]
    df_result = util.trimDF(df,
        delete_columns=[cn.PREDICTOR, cn.FREQ_PREDICTOR])
    return df_result


if __name__ == '__main__':
  dm = DataModel()
  dm.do()
