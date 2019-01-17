"""Data related constants."""

from microbepy.common.schema import TableSchemas

import numpy as np
import os
from pathlib import Path


# Database
DB_NAME = "microbepy"

# Species
SPECIES_MIX_DVH = 'D'
SPECIES_MIX_MMP = 'M'
SPECIES_MIX_BOTH = 'B'  # One M and one D
SPECIES_MIX_EPD = 'E'  # Endpoint dilution
SPECIES_MIX_LINE = 'L'
SPECIES_MIX_UNKNOWN = 'U'
SPECIES_MIXES = [SPECIES_MIX_DVH, SPECIES_MIX_MMP, 
    SPECIES_MIX_BOTH, SPECIES_MIX_EPD, SPECIES_MIX_LINE,
    SPECIES_MIX_UNKNOWN]
SPECIES_DICT = {
    SPECIES_MIX_DVH: "DVH",
    SPECIES_MIX_MMP: "MMP",
    }

# Transfers
TRANSFER_1000G = 152
TRANSFER_DEFAULT = TRANSFER_1000G

# Lines
LINE_WT = 'WT'
LINE_HA2 = 'HA2'
LINE_HR2 = 'HR2'
LINE_AN1 = 'AN1'  # AN_Dv-Ancestor-1, AN_Mm-Ancestor-1
LINE_AN2 = 'AN2'  # AN_Dv-Ancestor-2, AN_Mm-Ancestor-2
LINE_ANC = 'ANC'  # AN_Coculture-Ancestor
LINE_AN = 'AN'  # TODO: Is this needed
LINE_ALL = 'all'
LINE_UE3 = "UE3"
LINE_CIS = [LINE_HA2, LINE_HR2, LINE_UE3]

# Mutations
MUTATION_SEPARATOR = "."  # Separates components of a mutation key

# Isolates
ISOLATE_DEFAULT = '*'  # Default value for an isolate component
ISOLATE_SEPARATOR = "."  # Separates components of an isolate key
ISOLATE_UNKNOWN = "isolate_unknown"
#
ISOLATE_LINE = "isolate_line"
ISOLATE_EPD = "isolate_epd"  # EPD isolate
ISOLATE_CLONE = "isolate_clone"
ISOLATE_ANC = "isolate_anc"
ISOLATE_ANX = "isolate_anx"
ISOLATE_PAIR = "isolate_pair"  # (dvh, mmp)

# Experiments
EXPERIMENT_CI = "CI"
EXPERIMENT_SC = "SC"
EXPERIMENT_NULL = ""

# None values
NONE = np.nan

# Column in the Data Model
AA_CHANGE = 'aa_change'  # str amino acid
ACCESSION = "accession"  # str  reference number
CHANGED_NT = 'changed_nt'  # str nucleotide sequence
CLONE = "clone"  # Number of the clone
CODON_CHANGE = 'codon_change'  # str codon
COG = "cog"
COG_DESC = "cog_desc"
COG_FUN = "cog_fun"
COMMUNITY = "community"  # Value ISOLATE_{CLONE, EPD, LINE}
EC = "ec"
EC_DESC = "ec_desc"
EFFECT = 'effect'  # str HIGH, MODERATE, MODIFIER, LOW
EPD_ID = 'endpoint_dilution_id'  # int
EXPERIMENT = "experiment"  # str SC (single cell), CI (clonal isolate)
GATK = "gatk"  # float percent predicted by GATK
GENE_DESC = "gene_desc"  # phrase describing the GENE function
GENE_GI = "gene_gi"
GENE_ID = 'gene_id'  # str name of gene
GENE_NAME = "gene_name"
GENE_POSITION = "gene_position"  # <gene>.<position>
GGENE_ID = 'ggene_id'  # str name of gene if non-null; else gene_position
GO = "go"
INITIAL_NT = 'initial_nt'  # str original nucleotide sequence
INITIAL_CODON = 'initial_codon'  # str original codon
IS_LOW_COVERAGE_ISOLATE = 'is_low_coverage_isolate'  # bool
IS_AN_MUTATION = 'is_an_mutation'  # bool
LINE = 'line'
KEY_CULTURE = 'key_culture'
KEY_ISOLATE = 'key_isolate'  # str 6 component isolate specification
KEY_MUTATION = 'key_mutation'
LOCUS_ID = "locus_id"
MUTATION = 'mutation'
POSITION = 'position'
PREDICTOR = "predictor"
RATE = 'rate'
READ_NUMBER = "read_number"
RELATIVE_POSITION = 'relative_position'  # relative position of mutation within the gene
SAMTOOLS = "samtools"
SCAFFOLD_ID = "scaffold_id"
SOURCE = 'source'
SPECIES = 'species'
SPECIES_MIX = 'species_mix'
START = 'start'  # Start BP of a gene
STOP = 'stop'  # Stop position of a gene
STRAND = 'strand'  # Direction of the strand
TIGR_FAM = "tigr_fam"
TIGR_ROLES = "tigr_roles"
TRANSFER = "transfer"  # Transfer number
TYPE = 'type'
VARIANT_ID = "variant_id"
VARSCAN = "varscan"
YIELD = 'yield'

# Columns in files used as inputs to create data model
FREQ = "freq"
FREQ_PREDICTOR = "freq_predictor"
KEY_ISOLATE_IMPUTE = 'key_isolate_impute'  # imputed for Low coverage isolate
ISOLATE1 = "isolate_1"
ISOLATE2 = "isolate_2"
SAMPLE = "sample"

# Other column names commonly used
ACC = "acc"  # model accuracy
AVG = "avg"  # Average
CATEGORICAL = "categorical"
CENTER = "center"
COUNT = 'count'
COUNT1 = 'count1'
COUNT2 = 'count2'
COUNT_ISOLATES = "count_isolates"
COUNT_MUTATIONS = "count_mutations"
DEPVAR = "depvar"  # Dependent variable
DEPVARS = [RATE, YIELD]
ESTIMATE = "estimate"  # Estimate from regression
FOLD = "fold"
FRACTION = 'fraction'  # A fraction
GROUP = "group"
IDX = 'idx'  # Used to count instances
INDEX = 'index'  # Created by pandas
ISOLATES = 'isolates'  # list of isolates
KEY_ISOLATE_DVH = "key_isolate_dvh"
KEY_ISOLATE_MMP = "key_isolate_mmp"
MUTATIONS = "mutations"
MUTATION_COLUMN = "mutation_column"
MUTATION_LENGTH = 'mutation_length'
LEAF_NODE = "leaf_node"  # Leaf node in a tree
LOWER = "lower"
MAX = "max"
MIN = "min"
MUTATION_COLUMNS = [GENE_ID, GGENE_ID, POSITION,
    KEY_MUTATION]
NULL_SET = set([])
OBSERVED = "observed"  # observed values
PARAMETER = "parameter"
PREDICT = "predict"
PREDVAR = "predvar"  # predicting variable
PROBABILITY = "probability"
REPLICATION = "replication"  # replication number
RESIDUAL = "residual"  # Residual from regression
RESIDUALSTD = "residualstd"  # Residual from regression in units of standard deviation
RNG = "rng"  # Range of values
RSQ = "rsq"  # R2
SIZE = "size"
SORT = "sort"  # Used for sorting
SSQ = "ssq"  # sum of squares
SET = 'set'
SIGLVL = "siglvl"
SLOPE = "slope"
STD = "std"  # Standard deviation
TEST_X = "test_x"
TEST_Y = "test_y"
TOTAL = 'total'
TSTAT = "tstat"  # t statistic
TRAIN_X = "train_x"
TRAIN_Y = "train_y"
UPPER = "upper"
VALUE = "value"
VARIABLE = "variable"
XAXIS = 'xaxis'

# Columns in the raw mutation data
MUTATION_RAW_COLUMNS = [VARIANT_ID,
    SOURCE, POSITION, INITIAL_NT, CHANGED_NT, EFFECT,
    TYPE, AA_CHANGE, GENE_ID,
    PREDICTOR, FREQ_PREDICTOR, READ_NUMBER]

# Table schemas in the data model
TABLE_CULTURE = "culture"
TABLE_MUTATION = "mutation"
TABLE_ISOLATE = "isolate"
TABLE_CULTURE_ISOLATE_LINK = "culture_isolate_link"
TABLE_ISOLATE_MUTATION_LINK = "isolate_mutation_link"
TABLE_GENE_DESCRIPTION = "gene_description"
TABLE_CULTURE_ISOLATE_MUTATION = "culture_isolate_mutation"
TABLE_GENOTYPE = "genotype"
TABLE_GENOTYPE_PHENOTYPE = "genotype_phenotype"

# Schemas
TABLE_SCHEMAS = TableSchemas()
# Define the columns and their type
TABLE_SCHEMAS.column_schemas.addSchema(
    [IS_AN_MUTATION, IS_LOW_COVERAGE_ISOLATE], bool)
TABLE_SCHEMAS.column_schemas.addSchema(
    [RATE, YIELD, GATK, VARSCAN, SAMTOOLS, FREQ], float)
TABLE_SCHEMAS.column_schemas.addSchema(
    [POSITION, EPD_ID, CLONE, TRANSFER, START, STOP, ], int)
TABLE_SCHEMAS.column_schemas.addSchema(
    [KEY_CULTURE, SPECIES_MIX, KEY_ISOLATE, LINE, 
    SPECIES, EXPERIMENT, MUTATION, GENE_POSITION,
    LOCUS_ID, ACCESSION, GENE_GI, SCAFFOLD_ID, 
    STRAND, GENE_NAME, GENE_DESC, COG, COG_FUN, COG_DESC, TIGR_FAM,
    TIGR_ROLES, GO, EC, EC_DESC,
    SOURCE, INITIAL_NT, CHANGED_NT, EFFECT, TYPE, READ_NUMBER,
    AA_CHANGE, GENE_ID, GGENE_ID, KEY_MUTATION,
    ], str)
# Specify the columns in each table
TABLE_SCHEMAS.addSchema(TABLE_CULTURE,
    [KEY_CULTURE, RATE, YIELD, SPECIES_MIX],
    key=[KEY_CULTURE])
TABLE_SCHEMAS.addSchema(TABLE_ISOLATE,
    [KEY_ISOLATE, LINE, TRANSFER, EPD_ID, CLONE, SPECIES, EXPERIMENT, 
    IS_LOW_COVERAGE_ISOLATE],
    [KEY_ISOLATE])
TABLE_SCHEMAS.addSchema(TABLE_MUTATION,
    [SOURCE, POSITION, INITIAL_NT, CHANGED_NT, 
    EFFECT, TYPE, AA_CHANGE, GENE_ID, GGENE_ID,
    IS_AN_MUTATION, GENE_POSITION,
    MUTATION, KEY_MUTATION],
    [KEY_MUTATION])
TABLE_SCHEMAS.addSchema(TABLE_ISOLATE_MUTATION_LINK,
    [KEY_ISOLATE, KEY_MUTATION,
    GATK, VARSCAN, SAMTOOLS, READ_NUMBER, FREQ],
    [KEY_ISOLATE, KEY_MUTATION])
TABLE_SCHEMAS.addSchema(TABLE_CULTURE_ISOLATE_LINK,
    [KEY_ISOLATE, KEY_CULTURE],
    [KEY_ISOLATE, KEY_CULTURE])
TABLE_SCHEMAS.addSchema(TABLE_GENE_DESCRIPTION,
    [GENE_ID, LOCUS_ID, ACCESSION, GENE_GI, SCAFFOLD_ID, START, STOP,
    STRAND, GENE_NAME, GENE_DESC, COG, COG_FUN, COG_DESC, TIGR_FAM,
    TIGR_ROLES, GO, EC, EC_DESC],
    [])
columns = TABLE_SCHEMAS.getColumns([
    TABLE_CULTURE, TABLE_ISOLATE, TABLE_MUTATION,
    TABLE_CULTURE_ISOLATE_LINK, TABLE_ISOLATE_MUTATION_LINK,
    TABLE_GENE_DESCRIPTION])
TABLE_SCHEMAS.addSchema(TABLE_CULTURE_ISOLATE_MUTATION,
    columns, [KEY_CULTURE, KEY_ISOLATE, KEY_MUTATION])
columns = TABLE_SCHEMAS.getColumns([
    TABLE_MUTATION,
    TABLE_ISOLATE_MUTATION_LINK,
    TABLE_ISOLATE,
    TABLE_GENE_DESCRIPTION,
    ])
TABLE_SCHEMAS.addSchema(TABLE_GENOTYPE, 
    columns,
    [KEY_MUTATION, KEY_ISOLATE]
    )
columns = TABLE_SCHEMAS.getColumns([
    TABLE_MUTATION,
    TABLE_ISOLATE_MUTATION_LINK,
    TABLE_ISOLATE,
    TABLE_CULTURE_ISOLATE_LINK,
    TABLE_CULTURE,
    TABLE_GENE_DESCRIPTION,
    ])
TABLE_SCHEMAS.addSchema(TABLE_GENOTYPE_PHENOTYPE, 
    columns,
    [KEY_MUTATION, KEY_ISOLATE, KEY_CULTURE, TRANSFER]
    )
# Add functional dependencies
TABLE_SCHEMAS.addFD(KEY_MUTATION, GENE_POSITION)
TABLE_SCHEMAS.addFD(KEY_MUTATION, GGENE_ID)
TABLE_SCHEMAS.addFD(KEY_MUTATION, GENE_ID)
TABLE_SCHEMAS.addFD(GENE_POSITION, GENE_ID)
TABLE_SCHEMAS.addFD(GGENE_ID, GENE_ID)
TABLE_SCHEMAS.addFD(GENE_POSITION, GGENE_ID)

# Mutations
EFFECT_NONE = NONE
EFFECT_HIGH = 'HIGH'
EFFECT_MODERATE = 'MODERATE'
EFFECT_MODIFIER = 'MODIFIER'
EFFECT_LOW = 'LOW'
EFFECTS = [EFFECT_HIGH, EFFECT_MODERATE, 
    EFFECT_LOW, EFFECT_MODIFIER, EFFECT_NONE]

# CGI Matrices
CGI_VALUES_BINARY = 'cgi_values_binary'
CGI_VALUES_POSITION = 'cgi_values_position'

# Isolates
# The following isolates are excluded from the analysis because of
# low coverage for the reads
LOW_COVERAGE_ISOLATES = [
    'HA2.152.01.02.D.CI', 
    'HA2.152.01.03.D.CI', 
    'HA2.152.05.02.D.CI', 
    'HR2.152.10.03.D.CI', 
    'UE3.152.03.01.D.CI',
    'HA2.152.08.02.M.CI', 
    'HA2.152.08.03.M.CI', 
    'HA2.152.09.02.M.CI', 
    'UE3.152.03.02.M.CI',
    ]

# The gene_description data are incomplete and so
# several genes are not present that do have mutations
GENES_MISSING_IN_GENE_DESCRIPTION = ['DVU0109', 'DVU0436', 'DVU0845', 
    'DVU1214', 'DVU1571', 'DVU2456', 'DVU2955']

# Test related
IS_TEST = False
TEST_PATH = "/tmp"

# Plot parameters
PLT_FIGSIZE = "figsize"
PLT_LEGEND = "legend"
PLT_TITLE = "title"
PLT_XLABEL = "xlabel"
PLT_YLABEL = "ylabel"
PLT_XLIM = "xlim"
PLT_YLIM = "ylim"
PLT_XTICKLABELS = "xticklabels"
PLT_YTICKLABELS = "yticklabels"
PLT_XAXISTICKTOP = "xaxisticktop"
PLT_COLORBAR = "colorbar"

# Mutations in increasing order of refinement
MUTATION_COLUMNS = [GENE_ID, GGENE_ID, POSITION, KEY_MUTATION]

# Special constants
ONE = 1
ZERO = 0
BINARY_VALUES = [ZERO, ONE]

# Predictor transformation types
TRANSFORM_DEFAULT  \
    = "genome_model_transform_default"
TRANSFORM_NONE  \
    = "genome_model_transform_none"
TRANSFORM_LOW_FREQUENCY_ISOLATES  \
    = "genome_model_transform_low_frequency_isolates"
TRANSFORM_ONLY_LOW_FREQUENCY_ISOLATES  \
    = "genome_model_transform_only_low_frequency_isolates"

# Significance levels
SL_TSTAT = "sl_tstat"  # Significance level computed by t-statistic
SL_RESAMPLE = "sl_resample"  # Computed from resample

# Extrema
EXT_MIN = "min"
EXT_MAX = "max"
EXT_VALUES = [EXT_MIN, EXT_MAX]
      
# File paths
CONFIG_DIR = ".microbepy"  # Name of the directory with config info
HOME_DIR = str(Path.home())
CONFIG_DIR_PATH = os.path.join(HOME_DIR, CONFIG_DIR)
CONFIG_FILE = "config.yaml"  # Name of the configuration file
CONFIG_FILE_PATH = os.path.join(CONFIG_DIR_PATH, CONFIG_FILE)
SQLDB_PATH = None  # Default value for the SQLDB_PATH
SQLDB_PATH_NAME = "SQLDB_PATH"
YAML_DEFAULT = {SQLDB_PATH_NAME: SQLDB_PATH}
SQLDB_FILE = "%s.db" % DB_NAME
