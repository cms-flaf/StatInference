import copy
import math
import re
from enum import Enum

from ..common.tools import importROOT
ROOT = importROOT()

from CombineHarvester.CombineTools.ch import SystMap

class UncertaintyScale(Enum):
  Up = 1
  Down = -1

class UncertaintyType(Enum):
  auto = 0
  lnN = 1
  shape = 2

class Uncertainty:
  def __init__(self, name, processes=None, eras=None, channels=None, categories=None):
    self.name = name
    self.processes = processes
    self.eras = eras
    self.channels = channels
    self.categories = categories

  def appliesTo(self, process, era, channel, category):
    match_subprocesses = []
    if process.subprocesses:
      for subpr in process.subprocesses:
        match_subprocesses.append(Uncertainty.hasMatch(subpr, process.subprocesses))
    match_one_subpr=False
    for item in match_subprocesses:
      if item==True:
        match_one_subpr=True
        break
    #print(Uncertainty.hasMatch(process.name, self.processes) or match_one_subpr)
    #print(Uncertainty.hasMatch(era, self.eras))
    #print(Uncertainty.hasMatch(channel, self.channels))
    #print(Uncertainty.hasMatch(category, self.categories))
    return  (Uncertainty.hasMatch(process.name, self.processes) or match_one_subpr) \
        and Uncertainty.hasMatch(era, self.eras) \
        and Uncertainty.hasMatch(channel, self.channels) \
        and Uncertainty.hasMatch(category, self.categories)

  def resolveType(self, nominal_shape, shape_variations, auto_lnN_threshold, asym_value_threshold):
    return self

  @staticmethod
  def fitFlat(histogram):
    def constantFunction(x,par):
      return par[0]
    fit_func = ROOT.TF1("fit_func", constantFunction, 0, 10, 1)
    fit_func.SetParameter(0, 1.0)
    histogram.Fit(fit_func, "nq")
    chi2 = fit_func.GetChisquare()
    ndf = fit_func.GetNDF()
    p_value = ROOT.TMath.Prob(chi2, ndf)
    fit_param = fit_func.GetParameter(0)
    fit_param_error = fit_func.GetParError(0)
    return chi2, p_value, fit_param, fit_param_error

  @staticmethod
  def haveCompatibleShapes(hist_a, hist_b, p_thr):
    hist_b = hist_b.Clone()
    hist_b.Divide(hist_a)
    _, p_value, _, _ = Uncertainty.fitFlat(hist_b)
    return p_value > p_thr

  @staticmethod
  def hasMatch(value, patterns):
    if len(patterns) == 0:
      return True
    for pattern in patterns:
      #print(f"value = {value}")
      #print(pattern)
      #print(value == pattern)
      if pattern[0] == '^':
        if re.match(pattern, value):
          return True
      elif value == pattern:
        return True
    return False

  @staticmethod
  def fromConfig(entry):
    name = entry["name"]
    unc_type = UncertaintyType[entry["type"]]

    def getPatternList(key):
      if key not in entry:
        return []
      v = entry[key]
      if type(v) == str:
        return [ v ]
      return v

    args = {}
    for key in [ "processes", "eras", "channels", "categories"]:
      args[key] = getPatternList(key)

    if unc_type == UncertaintyType.lnN:
      value = entry["value"]
      if type(value) in [ float, int ]:
        value = float(value)
      elif type(value) == list and len(value) == 2:
        value = { UncertaintyScale.Down: value[0], UncertaintyScale.Up: value[1] }
      elif type(value) == dict :
        value = value
      else:
        raise RuntimeError("Invalid lnN uncertainty value")

    if unc_type == UncertaintyType.lnN:
      unc = LnNUncertainty(name, value, **args)
      unc.checkValue()
    elif unc_type == UncertaintyType.shape:
      unc = ShapeUncertainty(name, **args)
    elif unc_type == UncertaintyType.auto:
      unc = AutoUncertainty(name, **args)
    else:
      raise RuntimeError("Invalid uncertainty type")

    return unc

class LnNUncertainty(Uncertainty):
  def __init__(self, name, value, **kwargs):
    super().__init__(name, **kwargs)
    self.value = value

  @property
  def type(self):
    return UncertaintyType.lnN

  @property
  def needShapes(self):
    return False

  def canIgnore(self, value_thr):
    if type(self.value) == float:
      return abs(self.value) < value_thr
    return abs(self.value[UncertaintyScale.Up]) < value_thr and abs(self.value[UncertaintyScale.Down]) < value_thr

  def checkValue(self, raise_error=True):
    pass_zero_check = True
    pass_sign_check = True
    pass_magnitude_check = True
    if type(self.value) == float:
      if self.value == 0:
        pass_zero_check = False
      if abs(self.value) >= 1:
        pass_magnitude_check = False
    elif type(self.value) == dict:
      if self.value[UncertaintyScale.Down] == 0 or self.value[UncertaintyScale.Up] == 0:
        pass_zero_check = False
      if self.value[UncertaintyScale.Down] * self.value[UncertaintyScale.Up] > 0:
        pass_sign_check = False
      if abs(self.value[UncertaintyScale.Down]) >= 1 or abs(self.value[UncertaintyScale.Up]) >= 1:
        pass_magnitude_check = False
    else:
      raise RuntimeError("Invalid lnN uncertainty value")

    if not pass_zero_check and raise_error:
      raise RuntimeError(f"Value of lnN uncertainty {self.name} cannot be zero")
    if not pass_sign_check and raise_error:
      raise RuntimeError(f"Up/down values of lnN uncertainty {self.name} must have opposite signs")
    if not pass_magnitude_check and raise_error:
      raise RuntimeError(f"Value of lnN uncertainty {self.name} must be less than 100%")
    return pass_zero_check, pass_sign_check, pass_magnitude_check

  def valueToMap(self, digits=3):
    if type(self.value) == float:
      value = round(1 + self.value, digits)
    else:
      v_down = round(1 + self.value[UncertaintyScale.Down], digits)
      v_up = round(1 + self.value[UncertaintyScale.Up], digits)
      value = (v_down, v_up)
    return SystMap()(value)

class ShapeUncertainty(Uncertainty):
  @property
  def type(self):
    return UncertaintyType.shape

  @property
  def needShapes(self):
    return True

  def canIgnore(self, value_thr):
    return False

  def valueToMap(self, digits=3):
    return SystMap()(1.0)

class AutoUncertainty(Uncertainty):
  @property
  def needShapes(self):
    return True

  def canIgnore(self, value_thr):
    raise RuntimeError("Uncertainty type is not resolved yet")

  def valueToMap(self, digits=3):
    raise RuntimeError("Uncertainty type is not resolved yet")

  def resolveType(self, nominal_shape, shape_variations, auto_lnN_threshold, asym_value_threshold):
    nominal_shape = nominal_shape.Clone()
    for n in range(0, nominal_shape.GetNbinsX() + 2):
      nominal_shape.SetBinError(n, 0)
    all_compatible = True
    for unc_scale in [ UncertaintyScale.Up, UncertaintyScale.Down ]:
      if not Uncertainty.haveCompatibleShapes(nominal_shape, shape_variations[unc_scale], auto_lnN_threshold):
        all_compatible = False
        break
    args = { "processes": self.processes, "eras": self.eras, "channels": self.channels, "categories": self.categories }
    if all_compatible:
      unc_value = {}
      nominal_integral = nominal_shape.Integral()
      for unc_scale in [ UncertaintyScale.Up, UncertaintyScale.Down ]:
        scaled_integral = shape_variations[unc_scale].Integral()
        unc_value[unc_scale] = (scaled_integral - nominal_integral) / nominal_integral
      unc = LnNUncertainty(self.name, unc_value, **args)
      pass_zero_check, pass_sign_check, pass_magnitude_check = unc.checkValue(raise_error=False)
      if pass_sign_check and pass_magnitude_check:
        delta = abs(unc.value[UncertaintyScale.Up]) - abs(unc.value[UncertaintyScale.Down])
        if abs(delta) <= asym_value_threshold:
          max_value = max(abs(unc.value[UncertaintyScale.Up]), abs(unc.value[UncertaintyScale.Down]))
          unc.value = math.copysign(max_value, unc.value[UncertaintyScale.Up])
        return unc
    return ShapeUncertainty(self.name, **args)
