import copy
from ..common.param_parse import extractParameters, applyParameters

class Model:
  def getInputFileName(self, era, param_values=None):
    param_values = copy.deepcopy(param_values or {})
    param_values['ERA'] = era
    file_param_values = { }
    for param in self.file_parameters:
      file_param_values[param] = param_values[param]
    return applyParameters(self.input_file_pattern, file_param_values)

  def paramStr(self, param_values=None):
    # for nonresonant params: Ensure Model.paramStr(params) ignores 'mass'
    if param_values is None: return "nominal"
    clean = {k:v for k,v in param_values.items() if k != "mass"}
    if not clean: return "nominal"
    # for nonresonant params, e.g. kl1p0_kt1p0
    param_strs = [ f"{param}_{clean[param]}" for param in clean.keys() ]
    return '_'.join(param_strs)

  def massStr(self, param_values):
    # for mass: only called if 'mass' in params
    return str(int(param_values['mass']))   # standard CH mass axis

  @staticmethod
  def fromConfig(entry):
    model = Model()
    model.parameters = entry.get("parameters", [])
    model.param_dependent_bkg = entry.get("param_dependent_bkg", False)
    model.input_file_pattern = entry['input_file_pattern']
    model.file_parameters = extractParameters(model.input_file_pattern)
    file_model_pset = set(model.file_parameters) - { "ERA" }
    if len(file_model_pset - set(model.parameters)) > 0:
      raise RuntimeError("Invalid parameters in file name")
    if model.param_dependent_bkg:
      if len(model.parameters) == 0:
        raise RuntimeError("Parameter dependent background requires parameters")
      if len(model.parameters) != len(file_model_pset):
        raise RuntimeError("Parameter dependent background requires all parameters in the file name")
    return model
