""" 
Performance Studies for Data Acces modules.
Studies both memory and time.
"""

class PerfMakeDF(object):
  """ Performance studies for api.makeDF. """

  def __init__(self, csv_file):
    self.csv_file = csv_file  # Where results are written as CSV

  def analyze(self, num_table, num_column)
    """
    Reports on the time and storage for queries that access a specified
    number of tables with the num_columns per table.
    :param int num_table:
    :parm int num_column:
    """
