
import __init__
import util
from api import Api
import constants as cn
import api as aa
import helpers

import numpy as np
import os
import pandas as pd
import unittest
import random, string

API_OBJECT = Api()


def testConstructor():
  for _ in range(10):
    df = API_OBJECT.makeDF()
    
if __name__ == '__main__':
  testConstructor()
