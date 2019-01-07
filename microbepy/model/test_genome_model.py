import microbepy_init
import constants as cn
import helpers
import util
import genome_model as gm
from cv_regression import CVBinaryTreeRegression
from cv_classification import CVBinaryTreeClassification

from isolate import Isolate
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import unittest


IGNORE_TEST = False
    

class TestGenomeModel(unittest.TestCase):

  def setUp(self):
    self.cls = gm.GenomeModel
    self.flags = {'min_rsq': 0.05}

  def testConstructor(self):
    if IGNORE_TEST:
      return
    model = self.cls(cn.RATE, cn.GGENE_ID, **self.flags)
    self.assertIsNotNone(model.df_X)

  def makeANConstraint(self):
    return lambda r: (not Isolate.isAN(r[cn.KEY_ISOLATE_DVH]))  \
        and not (Isolate.isAN(r[cn.KEY_ISOLATE_MMP]))

  def testFit(self):
    if IGNORE_TEST:
      return
    def getLine(isolate_stgs):
      for stg in isolate_stgs:
        iso = Isolate.create(stg)
        if iso.line in cn.LINE_CIS:
          return iso.line
      raise ValueError("Could not find line for %s")
    #
    def makeDF(constraints=None):
      mute_column = cn.GGENE_ID
      model = self.cls(cn.RATE, mute_column,
          transform_type=gm.cn.TRANSFORM_NONE,
          cv_model_cls=CVBinaryTreeClassification,
          constraints=constraints)
      cvr = model.fit()
      self.assertTrue(helpers.isValidDataFrame(cvr.df_parameter,
          cvr.df_parameter.columns))
      self.assertTrue(helpers.isValidDataFrame(cvr.df_predict,
          cvr.df_predict.columns))
      return cvr
    def plotCvr(cvr, line):
      dff = cvr.df_predict
      plt.scatter(dff[cn.OBSERVED], dff[cn.ESTIMATE])
      plt.plot([-2, 2], [-2, 2], color="red")
      plt.xlim([-2, 2])
      plt.ylim([-2, 2])
      nparm = len(cvr.df_parameter) - 1
      if cn.RSQ in cvr.df_parameter:
        col = cn.RSQ
      else:
        col = cn.ACC
      plt.title("%s. SCORE1: %f, SCORE2: %f, nparm: %d" % 
          (line, cvr.df_parameter.loc[col, cn.AVG], 
          cvr.score, nparm))
      plt.show()
    #
    constraints=[self.makeANConstraint()]
    cvr = makeDF(constraints=constraints)
    #
    if False:
      print ("****All lines")
      print (cvr.score)
      print (cvr.df_parameter)
      print (cvr.fitted_models[0])
      import pdb; pdb.set_trace()
    #
    if False:
      for line in LINE_CIS:
        print ("****Line: %s" % line)
        constraints=[self.makeANConstraint(),
            lambda r: r[cn.LINE] == line]
        cvr = makeDF(constraints=constraints)
        plotCvr(cvr, line)
    if False:
      plotCvr(cvr, 'all')

  def testFitRegression(self):
    if IGNORE_TEST:
      return
    model = self.cls(cn.RATE, cn.GGENE_ID, num_folds=-1,
        cv_model_cls=CVBinaryTreeRegression)
    cvr = model.fit()
    count = cvr.df_parameter[cn.COUNT].tolist()[0]
    self.assertGreater(count, 30)

  def testFitClassifier(self):
    if IGNORE_TEST:
      return
    model = self.cls(cn.RATE, cn.GGENE_ID, num_folds=-1,
        cv_model_cls=CVBinaryTreeClassification)
    cvc = model.fit()
    trues = [int(v) in [0, 1] for v in model.df_y[model.col_y]]
    self.assertTrue(all(trues))
    count = cvc.df_parameter[cn.COUNT].tolist()[0]
    self.assertGreater(count, 30)

  def testRegressByLine(self):
    if IGNORE_TEST:
      return
    result = self.cls.doByLine(cn.LINE_CIS, cn.RATE, cn.GGENE_ID,
        min_incr_score=0)
    if False:
      for line in cn.LINE_CIS:
        print ("\n\n***%s" % line)
        print (result[line].score)
        print (result[line].df_parameter)
    self.assertEqual(set(cn.LINE_CIS), set(result.keys()))
    for line in cn.LINE_CIS:
      df = result[line].df_parameter
      self.assertTrue(isinstance(df, pd.DataFrame))

  def testMakeClassificationData(self):
    if IGNORE_TEST:
      return
    genome_model = self.cls(cn.RATE, cn.GGENE_ID,
        percentile_threshold=50)
    genome_model.fit()
    length = len(genome_model.df_X)
    #
    genome_model = self.cls(cn.RATE, cn.GGENE_ID,
        percentile_threshold=25)
    genome_model.fit()
    self.assertGreater(length, len(genome_model.df_X))
    trues = [int(v) in [0, 1] 
        for v in genome_model.df_y[genome_model.col_y]]
    self.assertTrue(all(trues))
    #
    genome_model = self.cls(cn.RATE, cn.GGENE_ID,
        cv_model_cls=CVBinaryTreeRegression,
        percentile_threshold=25)
    genome_model.fit()
    trues = [int(v) not in [0, 1] 
        for v in genome_model.df_y[genome_model.col_y]]
    self.assertTrue(any(trues))
    #


if __name__ == '__main__':
    unittest.main()
