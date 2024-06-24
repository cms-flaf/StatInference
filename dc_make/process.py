from ..common.param_parse import extractParameters, applyParameters, parameterListToDict

class Process:
  def __init__(self, name, hist_name, is_signal=False, is_data=False, is_asimov_data=False, scale=1, params=None, allowNegativeContributions=False,subprocesses=[], isRelatedTo=None):
    self.name = name
    self.hist_name = hist_name
    self.is_signal = is_signal
    self.is_data = is_data
    self.is_asimov_data = is_asimov_data
    self.is_background = not (is_signal or is_data)
    self.scale=scale
    self.params = params
    self.allowNegativeContributions = allowNegativeContributions
    self.isRelatedTo = isRelatedTo
    self.subprocesses = subprocesses

    if is_data and is_signal:
      raise RuntimeError("Data and signal flags cannot be set simultaneously")
    if is_asimov_data and not is_data:
      raise RuntimeError("Asimov data flag can only be set for data processes")

    if is_signal:
      self.type = "signal"
    elif is_data:
      self.type = "data"
    else:
      self.type = "background"

  def __str__(self):
    str_rep =  f"Process({self.name}, type={self.type}, subprocesses={self.subprocesses}"
    if self.is_signal:
      str_rep += f", params={self.params}"
    str_rep += ")"
    return str_rep

  def hasCompatibleModelParams(self, model_params, param_dependent_bkg):
    if self.is_data or self.is_background:
      if model_params and not param_dependent_bkg:
        return False
      return True
    if not model_params:
      return False
    for param, value in self.params.items():
      if model_params[param] != value:
        return False
    return True

  @staticmethod
  def fromConfig(entry, model):
    if type(entry) == str:
      return [ Process(entry, entry) ]
    if type(entry) != dict:
      raise RuntimeError("Invalid entry type")
    base_name = entry["process"]
    base_hist_name = entry.get("hist_name", base_name)
    is_signal = entry.get("is_signal", False)
    is_data = entry.get("is_data", False)
    is_asimov_data = entry.get("is_asimov_data", False)
    allowNegativeContributions = entry.get("allowNegativeContributions", False)
    isRelatedTo = entry.get("isRelatedTo", False)
    subprocesses = entry.get("subprocesses", [])
    scale = entry.get("scale", 1)
    if type(scale) == str:
      scale = eval(scale)
    if 'param_values' not in entry:
      if is_signal and len(model.parameters) > 0:
        raise RuntimeError("Signal process must have parameter values")
      return [ Process(base_name, base_hist_name, is_signal=is_signal, is_data=is_data, is_asimov_data=is_asimov_data,
                       scale=scale,allowNegativeContributions=allowNegativeContributions,subprocesses=subprocesses, isRelatedTo=isRelatedTo) ]

    parameters = model.parameters if is_signal else extractParameters(base_name)
    param_values = entry["param_values"]
    if type(param_values) != list or len(param_values) == 0:
      raise RuntimeError("Invalid parameter values")
    processes = []
    for param_entry in param_values:
      param_dict = parameterListToDict(parameters, param_entry)
      name = applyParameters(base_name, param_dict)
      hist_name = applyParameters(base_hist_name, param_dict)
      processes.append(Process(name, hist_name, is_signal=is_signal, is_data=is_data, is_asimov_data=is_asimov_data,
                               scale=scale, params=param_dict))
    return processes