"""
Utilities for manipulating isolates and "generalized" isolates (e.g., lines).

Key conepts:
  line (or ancestral line) - the experimental conditions under which an incubation is done (resistance + stiring)
  line replica - instances of the same experimental conditions
  transfer - the transfer number from which the isolate was obtained. There is a linear relationship
      between the transfer number and the generation with transfer 152 being 1,000 generations.
  endpoint dilution (EPD) - a sample taken from a line replica after 1K generations
  endpoint dilution ID - a numeric identifier for an EPD
  isolate - a single genotype
  cell - a single organism
  species - each isolate is either a DVH or an MMP
  clone - an integer identifier of an isolate for an endpoint dilution
  
Isolates are encoded as follows: LLR.TTT.PP.CC.S.EE
  LLR - line replica (so R is a single digit). 
  TTT - transfer number, an integer up to 3 digits
  PP - two digit endpoint dilution
  CC - 2 digit clone number
  S -  1 character species identifier ('D', 'M')
  EE - 2 character experiment. CI - clonal isolate; SC - single cell

An EPD community has the format: LLR.TTT.PP

An clone pairing ID has the format: LLR.TTT.PP.CC

Wild type isolates begin with the string 'WT' and ancestral types with 'AN' followed by 'S' with '*' in the other positions.
"""

import microbepy_init
import util
import constants as cn

import copy
import os
import numpy as np
import pandas as pd

NONWT_ISOLATE_LENGTH = 14
WT_ISOLATE_LENGTH = 6


##################### HELPER FUNCTIONS ####################
def checkDefaults(default_values, non_default_values):
  """
  Checks that conditions hold for a set of values.
  :param list-of-str default_values: 
      values that should be cn.ISOLATE_DEFAULT
  :param list-of-str non_default_values: 
      values that should not be cn.ISOLATE_DEFAULT
  """
  result = True
  defaults = set([cn.ISOLATE_DEFAULT, str(np.nan)])
  if not all([str(x) in defaults
              for x in default_values]):
    result = False
  if not all([not str(x) in defaults
              for x in non_default_values]):
    result = False
  return result


##############################################
# Isolate class
##############################################
class Isolate(object):

  # Description of the different types of isolates. Isolates
  # can be classified by which components are set to the
  # default values. This classification is
  # used to validate isolates when they are created and to
  # classify them.
  TYPE_DICT = {

      # Ex: HA2.152.01.01.D.CI or
      # Ex: HA2.152.01.01.D.SC
      cn.ISOLATE_CLONE: lambda s: checkDefaults(
      [],
      [s.line, s.transfer, s.epd_id, s.clone, s.species, 
      s.experiment]),

      # Ex: HA2.152.01.*.*.*
      # Ex: HA2.152.01.02.*.*
      cn.ISOLATE_EPD: lambda s: checkDefaults(
      [s.species],
      [s.line, s.transfer, s.epd_id]),

      # Ex: HA2.12.*.*.*.*
      # Ex: WT.*.*.01.D.*
      # Ex: AN.*.*.01.D.*
      cn.ISOLATE_LINE: lambda s: checkDefaults(
      [s.epd_id, s.clone, s.species, s.experiment],
      [s.line, s.transfer])  \
      or (checkDefaults(
      [s.epd_id, s.transfer, s.experiment],
      [s.line, s.clone, s.species]) 
      and (s.line == cn.LINE_WT))
      or checkDefaults(
      [s.epd_id, s.transfer, s.experiment],
      [s.line, s.clone, s.species]) 
      and (s.line == cn.LINE_AN),

      # Ancestral co-culture: ANC.*.*.*.*.*
      cn.ISOLATE_ANC: lambda s: checkDefaults(
      [s.epd_id, s.transfer, s.clone, s.species, s.experiment],
      [s.line])  \
      and (s.line == cn.LINE_ANC),

      # Ancestral co-culture: AN1.*.*.*.D.*
      cn.ISOLATE_ANX: lambda s: checkDefaults(
      [s.epd_id, s.transfer, s.clone, s.experiment], 
      [s.line, s.species])  \
      and (s.line in [cn.LINE_AN1, cn.LINE_AN2]),

      # Unknown isolate: *.*.*.*.*.*
      # The species and experiment may or may not be known.
      cn.ISOLATE_UNKNOWN: lambda s: checkDefaults(
      [s.line, s.epd_id, s.transfer, s.clone],
      []),
      }

  def __init__(self,
      line=cn.ISOLATE_DEFAULT, 
      transfer=cn.ISOLATE_DEFAULT, 
      epd_id=cn.ISOLATE_DEFAULT, 
      clone=cn.ISOLATE_DEFAULT, 
      species=cn.ISOLATE_DEFAULT, 
      experiment=cn.ISOLATE_DEFAULT):

    def validate(value, func, msg):
      if value == cn.ISOLATE_DEFAULT:
        return
      elif func(value):
        return
      else:
        raise ValueError(msg)

    # line
    validate(line, 
        lambda x: (len(x) in [2, 3, 4]),
        "%s is an invalid line" % line)
    self.line = line
    # transfer
    validate(transfer, 
        lambda x: isinstance(int(x), int),
        "%s is an invalid transfer" % transfer)
    self.transfer = transfer
    # epd_id
    validate(epd_id, 
        lambda x: isinstance(int(x), int),
        "%s is an invalid epd_id" % epd_id)
    self.epd_id = epd_id
    # clone
    validate(clone, 
        lambda x: isinstance(int(x), int),
        "%s is an invalid clone" % clone)
    self.clone = clone
    # species
    validate(species, 
        lambda x: x in [cn.SPECIES_MIX_DVH, cn.SPECIES_MIX_MMP],
        "%s is an invalid species" % species)
    self.species = species
    # experiment
    self.experiment=experiment
    # Validate have consistent settings
    self._validateDefaultValues()

  def _validateDefaultValues(self):
    """
    Verifies the consistency of the assignment of default values to instance
    variables.
    :raises ValueError:
    """
    cls = self.__class__
    count = 0
    for _,f in cls.TYPE_DICT.items():
      if f(self):
        count += 1
    if count == 0:
      raise ValueError("%s does not match any isolate type" % str(self))
    elif count == 1:
      return
    else:
      raise ValueError("%s matches multiple isolate types" % str(self))

  @classmethod
  def create(cls, isolate_string):
    """
    Constructs an isolate object from an isolate string.
    :param str isolate_string:
    :return Isolate:
    """
    if util.isNull(isolate_string):
      return Isolate(
          line=cn.ISOLATE_DEFAULT, 
          transfer=cn.ISOLATE_DEFAULT, 
          epd_id=cn.ISOLATE_DEFAULT,
          clone=cn.ISOLATE_DEFAULT, 
          species=cn.ISOLATE_DEFAULT, 
          experiment=cn.ISOLATE_DEFAULT)
    elements = isolate_string.split(cn.ISOLATE_SEPARATOR)
    line = elements[0]
    transfer = cn.ISOLATE_DEFAULT 
    epd_id = cn.ISOLATE_DEFAULT 
    clone = cn.ISOLATE_DEFAULT 
    species = cn.ISOLATE_DEFAULT 
    experiment = cn.ISOLATE_DEFAULT
    
    # E.G. WT.D01
    if len(elements) == 2:
      species_clone = elements[1]
      species = species_clone[0]
      clone = species_clone[1:]
    # E.G., HA1.152.02.D01
    elif len(elements) == 4:
      transfer = int(elements[1])
      epd_id = elements[2]
      species = elements[3][0]
      if elements[3] == cn.ISOLATE_DEFAULT:
        clone = cn.ISOLATE_DEFAULT
        experiment = cn.ISOLATE_DEFAULT
      else:
        clone = elements[3][1:]
        experiment = cn.EXPERIMENT_CI
    # E.g., HA1.152.02.01.D.CI
    #  or AN.*.*.*.D.*
    #  or WT.*.*.01.D.*
    elif len(elements) == 6:
      transfer = elements[1]
      epd_id = elements[2]
      clone = elements[3]
      species = elements[4]
      experiment = elements[5]
    else:
      raise ValueError("%s is an invalid isolate string" % isolate_string)
    return Isolate(line=line, transfer=transfer, epd_id=epd_id,
        clone=clone, species=species, experiment=experiment)

  def getCommunity(self):
    """
    Determines the community specified by the isolate
    :return str: key of TYPE_DICT
    """
    cls = self.__class__
    for key in list(cls.TYPE_DICT.keys()):
      if cls.TYPE_DICT[key](self):
        return key
    raise RuntimeError("%s doesn't match any of the TYPE_DICT keys"
        % str(self))

  def getEPDCommunity(self):
    """
    Extracts the EPD Community from the isolate.
    :return str:
    """
    return "%s.%s.%s" % (self.line, self.transfer, self.epd_id)

  def getClonePairingID(self):
    """
    Extracts the EPD Community from the isolate.
    :return str:
    """
    return "%s.%s.%s.%s" % (self.line, self.transfer, self.epd_id, self.clone)

  def __str__(self):
    """
    String representation of an isolate.
    """
    return "%s.%s.%s.%s.%s.%s" % (
        self.line, str(self.transfer), self.epd_id, self.clone,
        self.species, self.experiment)

  @classmethod
  def isSpecies(cls, isolate_stg, species):
    if util.isNull(isolate_stg):
      return False
    isolate = cls.create(isolate_stg)
    return isolate.species == species

  @classmethod  
  def isEPD(cls, isolate_stg):
    """
    :param str isolate:
    :return bool: True if wild type isolate
                  False if non-wild type isolate
    """
    if util.isNull(isolate_stg):
      return False
    isolate = cls.create(isolate_stg)
    return isolate.getCommunity() == cn.ISOLATE_EPD
  
  @classmethod  
  def isAN(cls, isolate_stg):
    """
    :param str isolate_stg:
    :return bool: True if wild type isolate
                  False if non-wild type isolate
    :raises ValueError: not a valid isolate
    """
    if util.isNull(isolate_stg):
      return False
    isolate = cls.create(isolate_stg)
    return (isolate.line in 
        [cn.LINE_ANC, cn.LINE_AN, cn.LINE_AN1, cn.LINE_AN2]
    )
  
  @classmethod  
  def isSpecies(cls, isolate_stg, species):
    """
    Determines if the isolate is a particular species
    :param str isolate_stg:
    :param str species:
    :return bool: 
    """
    if util.isNull(isolate_stg):
      return False
    isolate = cls.create(isolate_stg)
    return isolate.species == species
  
  @classmethod  
  def isLine(cls, isolate_stg, line):
    """
    Determines if the isolate is a particular ancestral line.
    Handles case of line and combined line and line-replica.
    :param str isolate:
    :param str line:
    :return bool: 
    """
    if util.isNull(isolate_stg):
      return False
    isolate = cls.create(isolate_stg)
    offset = len(line)
    return isolate.line[0:offset] == line
