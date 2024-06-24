import itertools
import math
import os
import numpy as np
import yaml

from CombineHarvester.CombineTools.ch import CombineHarvester

from StatInference.common.tools import hasNegativeBins, listToVector, rebinAndFill, importROOT,hasRelevantNegativeBins
from .process import Process
from .uncertainty import Uncertainty, UncertaintyType, UncertaintyScale
from .model import Model
ROOT = importROOT()

def RenormalizeHistogram(histogram, norm, include_overflows=True):
    integral = histogram.Integral(0, histogram.GetNbinsX()+1) if include_overflows else histogram.Integral()
    if integral!=0:
        histogram.Scale(norm / integral)

def th1ToNumpy(histogram, include_overflow=False):
    bin_range = (0, histogram.GetNbinsX() + 1) if include_overflow else (1, histogram.GetNbinsX())
    n_bins = bin_range[1] - bin_range[0] + 1
    bin_contents = np.zeros(n_bins)
    bin_errors = np.zeros(n_bins)
    bin_numbers = np.zeros(n_bins, dtype=int)
    bin_indices = {}
    for bin_idx, bin_number in enumerate(range(bin_range[0], bin_range[1] + 1)):
        bin_contents[bin_idx] = histogram.GetBinContent(bin_number)
        bin_errors[bin_idx] = histogram.GetBinError(bin_number)
        bin_numbers[bin_idx] = bin_number
        bin_indices[bin_number] = bin_idx
    return bin_indices, bin_numbers, bin_contents, bin_errors

class NegativeBinSolution:
    def __init__(self):
        self.accepted = True
        self.integral = None
        self.has_zero_integral = False
        self.has_negative_integral = False
        self.negative_bins = set()
        self.negative_bins_within_error = set()
        self.relevant_negative_bins = set()
        self.changed_bins = set()

def resolveNegativeBins(histogram, allow_zero_integral=False, allow_negative_integral=False,
                        allow_negative_bins_within_error=False, max_n_sigma_for_negative_bins=1,
                        relevant_bins=None, include_overflow=False):
    relevant_bins = relevant_bins or []
    bin_indices, bin_numbers, bin_contents, bin_errors = th1ToNumpy(histogram, include_overflow)
    n_bins = len(bin_numbers)

    solution = NegativeBinSolution()
    solution.integral = np.sum(bin_contents)
    if solution.integral == 0:
        solution.has_zero_integral = True
        solution.accepted = allow_zero_integral
        return solution
    if solution.integral < 0:
        solution.has_negative_integral = True
        solution.accepted = allow_negative_integral
        return solution

    donor_integral = 0.
    negative_integral = 0.

    for bin_idx in range(n_bins):
        if bin_contents[bin_idx] >= 0:
            if bin_numbers[bin_idx] not in relevant_bins:
                donor_integral += bin_contents[bin_idx]
            continue
        negative_integral += bin_contents[bin_idx]
        solution.negative_bins.add(bin_numbers[bin_idx])
        #print(bin_idx, bin_numbers[bin_idx])
        zero_within_error = bin_contents[bin_idx] + max_n_sigma_for_negative_bins * bin_errors[bin_idx] >= 0
        if zero_within_error:
            solution.negative_bins_within_error.add(bin_numbers[bin_idx])
        if not allow_negative_bins_within_error or not zero_within_error:
            solution.accepted = False
        if bin_numbers[bin_idx] in relevant_bins:
            solution.relevant_negative_bins.add(bin_numbers[bin_idx])
            solution.accepted = False

    if abs(negative_integral) > donor_integral:
        solution.accepted = False

    if not solution.accepted or len(solution.negative_bins) == 0:
        return solution

    for bin_number in solution.negative_bins:
        bin_idx = bin_indices[bin_number]
        donor_bin_indices = []
        for i in range(n_bins):
            if i != bin_idx and bin_numbers[i] not in relevant_bins and bin_contents[i] > 0:
                donor_bin_indices.append(i)
        donor_bin_indices = sorted(donor_bin_indices, key=lambda i: abs(i-bin_idx))
        to_balance = -bin_contents[bin_idx]
        bin_errors[bin_idx] = math.sqrt(bin_contents[bin_idx] ** 2 + bin_errors[bin_idx] ** 2)
        bin_contents[bin_idx] = 0
        solution.changed_bins.add(bin_number)
        for donor_bin_idx in donor_bin_indices:
            donor_contribution = min(to_balance, bin_contents[donor_bin_idx])
            if donor_contribution <= 0: continue
            bin_errors[donor_bin_idx] = math.sqrt(donor_contribution ** 2 + bin_errors[donor_bin_idx] ** 2)
            bin_contents[donor_bin_idx] -= donor_contribution
            solution.changed_bins.add(bin_numbers[donor_bin_idx])
            to_balance -= donor_contribution
            if to_balance <= 0: break
        if to_balance > 0: # this should not happen, because we checked that negative_integral <= donor_integral
            raise RuntimeError(f"resolveNegativeBins: failed to balance bin {bin_number}")
    for bin_number in solution.changed_bins:
        bin_idx = bin_indices[bin_number]
        histogram.SetBinContent(bin_number, bin_contents[bin_idx])
        histogram.SetBinError(bin_number, bin_errors[bin_idx])

    return solution

class DatacardMaker:
  def __init__(self, cfg_file, input_path, hist_bins=None):
    self.cb = CombineHarvester()

    self.input_path = input_path
    with open(cfg_file, 'r') as f:
      cfg = yaml.safe_load(f)

    self.analysis = cfg["analysis"]
    self.eras = cfg["eras"]
    self.channels = cfg["channels"]
    self.categories = cfg["categories"]

    self.bins = []
    for era, channel, cat in self.ECC():
      bin = self.getBin(era, channel, cat, return_index=False)
      self.bins.append(bin)

    self.model = Model.fromConfig(cfg["model"])

    self.param_bins = {}
    self.processes = {}
    data_process = None
    has_signal = False
    for process in cfg["processes"]:
      new_processes = Process.fromConfig(process, self.model)
      for process in new_processes:
        if process.name in self.processes:
          raise RuntimeError(f"Process name {process.name} already exists")
        print(f"Adding {process}")
        self.processes[process.name] = process
        #self.processes[process.name]['subprocess'] = [process.subprocesses]
        if process.is_data:
          if data_process is not None:
            raise RuntimeError("Multiple data processes defined")
          data_process = process
        if process.is_signal:
          has_signal = True
          param_bin = self.model.paramStr(process.params)
          if param_bin in self.param_bins:
            raise RuntimeError(f"Signal process with parameters {param_bin} already exists")
          self.param_bins[param_bin] = process.params
    if data_process is None:
      raise RuntimeError("No data process defined")
    if not has_signal:
      raise RuntimeError("No signal process defined")

    self.uncertainties = {}
    for unc_entry in cfg["uncertainties"]:
      unc = Uncertainty.fromConfig(unc_entry)
      if unc.name in self.uncertainties:
        raise RuntimeError(f"Uncertainty {unc.name} already exists")
      self.uncertainties[unc.name] = unc

    self.autolnNThr = cfg.get("autolnNThr", 0.05)
    self.asymlnNThr = cfg.get("asymlnNThr", 0.001)
    self.ignorelnNThr = cfg.get("ignorelnNThr", 0.001)

    self.autoMCStats = cfg.get("autoMCStats", { 'apply': False })

    self.hist_bins = hist_bins
    if self.hist_bins is None:
      self.hist_bins = cfg.get("hist_bins", None)
    if self.hist_bins is not None:
      self.hist_bins = listToVector(self.hist_bins, 'double')

    self.input_files = {}
    self.shapes = {}

  def getBin(self, era, channel, category, return_name=True, return_index=True):
    name = f'{era}_{self.analysis}_{channel}_{category}'
    if not return_name and not return_index:
      raise RuntimeError("Invalid argument combination")
    if not return_index:
      return name
    index = self.bins.index(name)
    if not return_name:
      return index
    return (index, name)

  def cbCopy(self, param_str, process, era, channel, category):
    bin_idx, bin_name = self.getBin(era, channel, category)
    return self.cb.cp().mass([param_str]).process([process]).bin([bin_name])

  def ECC(self):
    return itertools.product(self.eras, self.channels, self.categories)

  def PPECC(self):
    param_bins = list(self.param_bins.keys())
    if not self.model.param_dependent_bkg:
      param_bins.append("*")
    return itertools.product(self.processes.keys(), param_bins, self.eras, self.channels, self.categories)

  def getInputFile(self, era, model_params):
    file_name = self.model.getInputFileName(era, model_params)
    if file_name not in self.input_files:
      full_file_name = os.path.join(self.input_path, file_name)
      file = ROOT.TFile.Open(full_file_name, "READ")
      if file == None:
        raise RuntimeError(f"Cannot open file {full_file_name}")
      self.input_files[file_name] = file
    return self.input_files[file_name]

  def prepareHist(self,file,hist_name,process):
    hist = file.Get(hist_name)
    if hist == None:
      raise RuntimeError(f"Cannot find histogram {hist_name} in {file.GetName()}")
    if self.hist_bins is not None:
      new_hist = ROOT.TH1F(hist.GetName(), hist.GetTitle(), len(self.hist_bins) - 1, self.hist_bins.data())
      rebinAndFill(new_hist, hist)
      hist = new_hist
    hist.SetDirectory(0)
    if process.scale != 1:
      hist.Scale(process.scale)
    return hist

  def getShape(self, process, era, channel, category, model_params,unc_name=None, unc_scale=None):
    key = (process.name, era, channel, category, unc_name, unc_scale)
    print(f"getting shape for {unc_name}")
    #print(f"process name {process.name}")
    if process.name == 'QCD' and category=='boosted':
      process.allowNegativeContributions = True
    if category == 'boosted': self.hist_bins = listToVector([245,600,3500], 'double')
    if key not in self.shapes:
      if process.is_data and (unc_name is not None or unc_scale is not None):
        raise RuntimeError("Cannot apply uncertainty to the data process")
      if process.is_asimov_data:
        hist = None
        for bkg_proc in self.processes.values():
          if bkg_proc.is_background:
            bkg_hist = self.getShape(bkg_proc, era, channel, category, model_params)
            if hist is None:
              hist = bkg_hist.Clone()
            else:
              hist.Add(bkg_hist)
        if hist is None:
          raise RuntimeError("Cannot create asimov data histogram")
      else:
        file = self.getInputFile(era, model_params)
        hist_name = f"{channel}/{category}/{process.hist_name}"
        hists = []
        if process.subprocesses:
          for subp in process.subprocesses:
            hist_name = f"{channel}/{category}/{subp}"
            if unc_name and unc_scale:
              hist_name += f"_{unc_name}{unc_scale}"
            hists.append(self.prepareHist(file,hist_name,process))
        else:
          hists.append(self.prepareHist(file,hist_name,process))
        hist = hists[0]
        if len(hists)>1:
          objsToMerge = ROOT.TList()
          for histy in hists:
            objsToMerge.Add(histy)
          hist.Merge(objsToMerge)
        hist.SetName(process.name)
        hist.SetTitle(process.name)
        relevant_bins = self.getRelevantBins(era, channel, category,unc_name, unc_scale=None)

        solution = resolveNegativeBins(hist,relevant_bins=relevant_bins, allow_zero_integral=category=='boosted')
        if not solution.accepted and process.allowNegativeContributions == False:
          axis = hist.GetXaxis()
          bins_edges = [ str(axis.GetBinLowEdge(n)) for n in range(1, axis.GetNbins() + 2)]
          bin_values = [ str(hist.GetBinContent(n)) for n in range(1, axis.GetNbins() + 1)]
          print(f'bins_edges: [ {", ".join(bins_edges)} ]')
          print(f'bin_values: [ {", ".join(bin_values)} ]')
          raise RuntimeError(f"Negative bins found in histogram {hist_name}")

      self.shapes[key] = hist
    return self.shapes[key]

  def getRelevantBins(self, era, channel, category,unc_name=None, unc_scale=None):
    relevant_bins = []
    for signal_proc in self.processes.values():
      if signal_proc.is_signal:
        model_params = signal_proc.params
        #param_str = self.model.paramStr(model_params)
        file = self.getInputFile(era, model_params)
        hist_name = f"{channel}/{category}/{signal_proc.hist_name}"
        hist = self.prepareHist(file,hist_name,signal_proc)
        axis = hist.GetXaxis()
        hist_integral = hist.Integral(1,axis.GetNbins() + 1)
        #print(axis.GetNbins())
        for nbin in range(1,axis.GetNbins() + 1):
          isImportant=False
          if hist_integral !=0 and hist.GetBinContent(nbin) / hist_integral >= 0.2:
            isImportant = True
          if len(relevant_bins) < axis.GetNbins():
            relevant_bins.append(isImportant)
          else:
            relevant_bins[nbin-1] = isImportant
    return relevant_bins

  def addProcess(self, proc, era, channel, category):
    bin_idx, bin_name = self.getBin(era, channel, category)
    process = self.processes[proc]

    def add(model_params, param_str):
      if process.is_data:
        self.cb.AddObservations([param_str], [self.analysis], [era], [channel], [(bin_idx, bin_name)])
      else:
        self.cb.AddProcesses([param_str], [self.analysis], [era], [channel], [proc], [(bin_idx, bin_name)], process.is_signal)
      shape = self.getShape(process, era, channel, category, model_params)
      shape_set = False
      def setShape(p):
        nonlocal shape_set
        print(f"Setting shape for {p}")
        if shape_set:
          raise RuntimeError("Shape already set")
        p.set_shape(shape, True)
        shape_set = True
      cb_copy = self.cbCopy(param_str, proc, era, channel, category)
      if process.is_data:
        cb_copy.ForEachObs(setShape)
      else:
        cb_copy.ForEachProc(setShape)

    if process.is_signal:
      model_params = process.params
      param_str = self.model.paramStr(model_params)
      add(model_params, param_str)
    elif self.model.param_dependent_bkg:
      for signal_proc in self.processes.values():
        if signal_proc.is_signal:
          model_params = signal_proc.params
          param_str = self.model.paramStr(model_params)
          add(model_params, param_str)
    else:
      add(None, "*")

  def addUncertainty(self, unc_name):
    unc = self.uncertainties[unc_name]
    for proc, param_str, era, channel, category in self.PPECC():
      process = self.processes[proc]
      if process.is_data: continue
      if not unc.appliesTo(process, era, channel, category): continue
      model_params = self.param_bins.get(param_str, None)
      if not process.hasCompatibleModelParams(model_params, self.model.param_dependent_bkg): continue

      nominal_shape = None
      shapes = {}
      #print(f"need shape? {unc.needShapes}")
      if unc.needShapes:
        model_params = self.param_bins.get(param_str, None)
        nominal_shape = self.getShape(self.processes[proc], era, channel, category, model_params)
        for unc_scale in [ UncertaintyScale.Up, UncertaintyScale.Down ]:
          shapes[unc_scale] = self.getShape(self.processes[proc], era, channel, category, model_params,
                                            unc_name, unc_scale.name)
      unc_to_apply = unc.resolveType(nominal_shape, shapes, self.autolnNThr, self.asymlnNThr)
      if unc_to_apply.canIgnore(self.ignorelnNThr):
        print(f"Ignoring uncertainty {unc_name} for {proc} in {era} {channel} {category}")
        continue
      systMap = unc_to_apply.valueToMap()
      cb_copy = self.cbCopy(param_str, proc, era, channel, category)
      cb_copy.AddSyst(self.cb, unc_name, unc_to_apply.type.name, systMap)
      if unc_to_apply.type == UncertaintyType.shape:
        shape_set = False
        def setShape(syst):
          nonlocal shape_set
          print(f"Setting unc shape for {syst}")
          if shape_set:
            raise RuntimeError("Shape already set")
          syst.set_shapes(shapes[UncertaintyScale.Up], shapes[UncertaintyScale.Down], nominal_shape)
          shape_set = True
        cb_copy = self.cbCopy(param_str, proc, era, channel, category).syst_name([unc_name])
        cb_copy.ForEachSyst(setShape)

  def writeDatacards(self, output):
    os.makedirs(output, exist_ok=True)
    background_names = [ proc_name for proc_name, proc in self.processes.items() if proc.is_background ]
    for proc_name, process in self.processes.items():
      if not process.is_signal: continue
      processes = [proc_name] + background_names
      param_list = [ self.model.paramStr(process.params) ]
      if not self.model.param_dependent_bkg:
        param_list.append('*')
      dc_file = os.path.join(output, f"datacard_{proc_name}.txt")
      shape_file = os.path.join(output, f"{proc_name}.root")
      self.cb.cp().mass(param_list).process(processes).WriteDatacard(dc_file, shape_file)

  def createDatacards(self, output, verbose=1):
    try:
      for era, channel, category in self.ECC():
        for process_name in self.processes.keys():
          self.addProcess(process_name, era, channel, category)
      for unc_name in self.uncertainties.keys():
        print(f"adding uncertainty: {unc_name}")
        self.addUncertainty(unc_name)
      if self.autoMCStats["apply"]:
        self.cb.SetAutoMCStats(self.cb, self.autoMCStats["threshold"], self.autoMCStats["apply_to_signal"],
                               self.autoMCStats["mode"])
      if verbose > 0:
        self.cb.PrintAll()
      self.writeDatacards(output)
    finally:
      for file in self.input_files.values():
        file.Close()