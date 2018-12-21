
import constants as cn

def readDataModelCSV(filename):
  """
  Reads the CSV for the file or table. Does fixups on data.
  :param str or TableSchema filename: filename with or without .csv extension
  :return pd.DataFrame:
  """
  def makeCSVFilename(ffile):
    name = "%s.csv" % ffile
    return name
  #
  csv_file = str(filename)
  if isinstance(filename, schema.TableSchema):
    csv_file = makeCSVFilename(filename.name)
  elif not filename.count(".csv") > 0:
    csv_file = makeCSVFilename(filename)
  path = getDataModelPath(csv_file)
  df = unifyNullValues(pd.read_csv(path))
  if cn.INDEX in df.columns:
    del df[cn.INDEX]
  pruneNullRows(df)
  return df

def makeCSVFilename(identifier):
  """
  :param TableSchema/str identifier:
  :return str: CSV file for the schema
  """
  if isStr(identifier):
    filename = identifier
  elif isinstance(identifier, schema.TableSchema):
    filename = identifier.name
  return "%s.csv" % filename

def makeDataframeFromXlsx(path, sheet_no=0):
  """
  Creates a dataframe from an xlsx file in a data directory
  :param str path: Path to data files
  :param int sheet_no: Number of the worksheet
  """
  data = pd.ExcelFile(path)
  return data.parse(sheet_no)
