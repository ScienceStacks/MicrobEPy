"""
Creates a CSV from an xlsx file.
If the xlsx has more than one sheet, then just creates a
CSV for the first sheet.
"""

import argparse
import pandas as pd

def createDataframe(in_path, out_path, worksheet=None):
  """
  Creates a CSV from the excel file
  :param str in_path: path to the excel file
  :param str out_path: path to the CSV file to create
  :param str worksheet: name of the worksheet in the excel file
  """
  data = pd.ExcelFile(in_path)
  if worksheet is None:
    worksheet = data.sheet_names[0]
  df = data.parse(worksheet=worksheet)
  df.to_csv(out_path, sep=",", index=False)

def extractArguments():
  """
  Extract the arguements
  :return dict: dictionary of argument, value pairs
  """
  parser = argparse.ArgumentParser(
      description='Convert xlsx to csv.')
  parser.add_argument('xlsx_file',
                      metavar='xlsx_file', 
                      type=str, nargs='+',
                      help='xlsx file')
  args = parser.parse_args()
  split_file = args.xlsx_file[0].split('.')
  last = len(split_file) - 1
  if split_file[last].lower() != 'xlsx':
    raise ValueError('Must be file type xlsx')
  filename = '.'.join(split_file[0:last])
  return {'filename': filename}

if __name__ == '__main__':
  arg_dict = extractArguments() 
  filename = arg_dict['filename']
  print("Processing file %s.xlsx..." % filename)
  in_file = "%s.xlsx" % filename
  out_file = "%s.csv" % filename
  createDataframe(in_file, out_file)
  print("  Created %s.csv" % filename)
