import itertools
import math
import os
import re
import yaml

from CombineHarvester.CombineTools.ch import CombineHarvester

from StatInference.common.tools import listToVector, rebinAndFill, importROOT, resolveNegativeBins, getRelevantBins
from .process import Process
from .uncertainty import Uncertainty, UncertaintyType, UncertaintyScale, MultiValueLnNUncertainty
from .model import Model
from .binner import Binner
ROOT = importROOT()

class DatacardMaker:
  def __init__(self, cfg_file, input_path, hist_bins=None, param_values=None):
    self.cb = CombineHarvester()

    self.input_path = input_path
    with open(cfg_file, 'r') as f:
      cfg = yaml.safe_load(f)

    self.analysis = cfg["analysis"]
    self.eras = cfg["eras"]
    self.channels = cfg["channels"]
    self.categories = cfg["categories"]
    self.signalFractionForRelevantBins = cfg['signalFractionForRelevantBins']

    self.bins = []
    for era, channel, cat in self.ECC():
      bin = self.getBin(era, channel, cat, return_index=False)
      self.bins.append(bin)

    self.model = Model.fromConfig(cfg["model"])
    #Define is_resonant: True only if there is exactly one parameter and it is "mass"
    ps = getattr(self.model, "parameters", [])
    self.is_resonant = (len(ps) == 1 and ps[0] == "mass")
    self.param_bins = {}
    self.processes = {}
    #Store per-(mass, process) parameters so it can be looked up correctly
    self.param_of = {}  # key: (mass_str, process_name) -> full params dict
    self.base_of = {}   # actual_proc_name -> base_process_name
    data_process = None
    has_signal = False
    for process in cfg["processes"]:
      if (type(process) != str) and process.get('is_signal', False):
        if param_values is not None:
          #previously "Overwriting signal masses" only for resonant case
          print(f"Overwriting signal parameters to {param_values}")
          process['param_values'] = param_values
      new_processes = Process.fromConfig(process, self.model)
      for process in new_processes:
        if process.name in self.processes:
          raise RuntimeError(f"Process name {process.name} already exists")
        print(f"Adding {process}")
        self.processes[process.name] = process
        if process.is_data:
          if data_process is not None:
            raise RuntimeError("Multiple data processes defined")
          data_process = process
        if process.is_signal:
          has_signal = True
          #?mass->nominal
          param_bin = self.model.paramStr(process.params)
          #Allow multiple signals to share same parameter point.
          self.param_bins.setdefault(param_bin, process.params)
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


    hist_bins = hist_bins or cfg.get("hist_bins", None)
    self.hist_binner = Binner(hist_bins)

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

  def cbCopy(self, mass_str, process, era, channel, category):
    bin_idx, bin_name = self.getBin(era, channel, category)
    return self.cb.cp().mass([mass_str]).process([process]).bin([bin_name])

  def ECC(self):
    return itertools.product(self.eras, self.channels, self.categories)

  def PPECC(self):
    param_bins = list(self.param_bins.keys())
    if not self.model.param_dependent_bkg:
      param_bins.append("*")
    return itertools.product(self.processes.keys(), param_bins, self.eras, self.channels, self.categories)

  def getInputFile(self, era, model_params, unc_name=None, unc_scale=None):
    file_name = self.model.getInputFileName(era, model_params)
    
    #Handle uncertainty files
    if unc_name and unc_scale:
      unc_file_name = file_name.replace("_Central.root", f"_{unc_name}.root")
      if os.path.exists(os.path.join(self.input_path, unc_file_name)):
        file_name = unc_file_name
        #print(f"DEBUG: Using uncertainty file: {file_name}")
      else:
        print(f"DEBUG: Uncertainty file not found, using Central file: {file_name}")
    
    if file_name not in self.input_files:
      full_file_name = os.path.join(self.input_path, file_name)
      file = ROOT.TFile.Open(full_file_name, "READ")
      if file == None:
        raise RuntimeError(f"Cannot open file {full_file_name}")
      self.input_files[file_name] = file
    return file_name, self.input_files[file_name]


  def getMultiValueLnUnc(self,unc,unc_name, process, era, channel, category, model_params):#, unc_name=None, unc_scale=None)
    file_name, file = self.getInputFile(era, model_params)
    hist_name = f"{channel}/{category}/{process.hist_name}"
    if unc.getUncertaintyForProcess(process.name) != None:
      return unc.getUncertaintyForProcess(process.name)
    elif process.subprocesses:
      unc_value_tot_down = 0.
      unc_value_tot_up = 0.
      yield_value_tot = 0.
      for subp in process.subprocesses:
        hist_name = f"{channel}/{category}/{subp}"
        subhist = file.Get(hist_name)
        #newhist = self.hist_binner.applyBinning(era, channel, category, model_params, subhist)
        if subhist == None:
          raise RuntimeError(f"Cannot find histogram {hist_name} in {file.GetName()}")
        axis = subhist.GetXaxis()
        yield_subproc = subhist.Integral(1,axis.GetNbins() + 1)
        unc_value = unc.getUncertaintyForProcess(subp)
        if unc_value != None:
          if yield_subproc == 0 : continue
          # print(unc_value)
          if isinstance(unc_value, dict):
            unc_value_tot_up += unc_value[UncertaintyScale.Up]*yield_subproc
            unc_value_tot_down += unc_value[UncertaintyScale.Down]*yield_subproc
          else:
            unc_value_tot_up += unc_value*yield_subproc
            unc_value_tot_down -= unc_value*yield_subproc
          yield_value_tot+=yield_subproc
      if unc_value_tot_up != 0. and unc_value_tot_down !=0 :
        return {UncertaintyScale.Down: unc_value_tot_down/yield_value_tot, UncertaintyScale.Up: unc_value_tot_up/yield_value_tot}
      return None
    return None


  def getShape(self, process, era, channel, category, model_params, unc_name=None, unc_scale=None):
    file_name, file = self.getInputFile(era, model_params, unc_name, unc_scale)
    #Add parameter tag to avoid cache key collisions across parameter points, eg kl_2_k2v_1
    param_tag = None if model_params is None else self.model.paramStr(model_params)
    key = (file_name, process.name, era, channel, category, unc_name, unc_scale, param_tag)
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
        base = f"{channel}/{category}/"
        
        hists = []
        if process.subprocesses:
          for subp in process.subprocesses:
            name = base + subp
            if unc_name and unc_scale:
              name += f"_{unc_name}_{unc_scale}"
            subhist = file.Get(name)
            if subhist is None:
              raise RuntimeError(f"Cannot find histogram {name} in {file.GetName()}")
            hists.append(self.hist_binner.applyBinning(era, channel, category, model_params, subhist))
        else:
          name = base + process.hist_name
          if unc_name and unc_scale:
            name += f"_{unc_name}_{unc_scale}"
          #print(f"DEBUG: Looking for histogram '{name}' in file '{file.GetName()}'")
          hist = file.Get(name)
          #print(f"DEBUG: Histogram result: {hist}")
          if hist is None:
            #file.ls()
            raise RuntimeError(f"Cannot find histogram {name} in {file.GetName()}")
          #print(f"DEBUG: Histogram found, applying binning...")
          hists.append(self.hist_binner.applyBinning(era, channel, category, model_params, hist))
        if len(hists) == 0:
          raise RuntimeError(f"hist list is empty for file {file.GetName()}")
        hist = hists[0]
        if len(hists)>1:
          for histy in hists[1:]:
            hist.Add(histy)
        hist.SetName(process.name)
        hist.SetTitle(process.name)

        hist.SetDirectory(0)
        if process.scale != 1:
          hist.Scale(process.scale)
        if process.is_signal:
            #store signal hists in self.shapes under the following tuple key:
            nominal_signal_key = ("signals", era, channel, category, param_tag)
            self.shapes.setdefault(nominal_signal_key, []).append(hist) 
        else:
          nominal_signal_key = ("signals", era, channel, category, param_tag)
          signal_processes_histograms = self.shapes.get(nominal_signal_key, [])
          relevant_bins = getRelevantBins(era, channel, category,signal_processes_histograms,self.signalFractionForRelevantBins,unc_name, unc_scale, model_params)
          solution = resolveNegativeBins(hist,relevant_bins=relevant_bins, allow_zero_integral=process.allow_zero_integral, allow_negative_bins_within_error=process.allow_negative_bins_within_error, max_n_sigma_for_negative_bins=process.max_n_sigma_for_negative_bins, allow_negative_integral=process.allow_negative_integral)

          if not solution.accepted:
            axis = hist.GetXaxis()
            bins_edges = [ str(axis.GetBinLowEdge(n)) for n in range(1, axis.GetNbins() + 2)]
            bin_values = [ str(hist.GetBinContent(n)) for n in range(1, axis.GetNbins() + 1)]
            bin_errors = [ str(hist.GetBinError(n)) for n in range(1, axis.GetNbins() + 1)]
            print(f'bins_edges: [ {", ".join(bins_edges)} ]')
            print(f'bin_values: [ {", ".join(bin_values)} ]')
            print(f'bin_errors: [ {", ".join(bin_errors)} ]')
            raise RuntimeError(
                f"Negative bins found in histogram for {channel}/{category}/{process.name}"
                + (f" (syst {unc_name}{unc_scale})" if unc_name and unc_scale else "")
            )
      self.shapes[key] = hist
    return self.shapes[key]




  def addProcess(self, proc, era, channel, category):
    bin_idx, bin_name = self.getBin(era, channel, category)
    process = self.processes[proc]
    #Modified add function to accept process_name parameter for unique signal names
    def add(model_params, mass_str, process_name):
      if process.is_data:
        self.cb.AddObservations([mass_str], [self.analysis], [era], [channel], [(bin_idx, bin_name)])
      else:
        self.cb.AddProcesses([mass_str], [self.analysis], [era], [channel], [process_name], [(bin_idx, bin_name)], process.is_signal)

      shape = self.getShape(process, era, channel, category, model_params)
      shape_set = False
      def setShape(p):
        nonlocal shape_set
        print(f"Setting shape for {p}")
        if shape_set:
          raise RuntimeError("Shape already set")
        p.set_shape(shape, True)
        shape_set = True
      cb_copy = self.cbCopy(mass_str, process_name, era, channel, category)
      if process.is_data:
        cb_copy.ForEachObs(setShape)
      else:
        cb_copy.ForEachProc(setShape)

    # for signals: iterate over each parameter point (mass, eft_tag)
    if process.is_signal:
      params = process.params
      if self.is_resonant:
        # Resonant: .mass = "<mass>", process name stays clean (e.g. XToHH)
        mass_str = str(int(params['mass']))
        actual_proc_name = process.name
      else:
        # Non-resonant: .mass = "*" and append EFT tag to the process name (e.g. ggHH_kl1p0_kt1p0)
        eft_tag = self.model.paramStr(params)   # must ignore 'mass'
        mass_str = "*"
        actual_proc_name = f"{process.name}_{eft_tag}"

      # Remember this parameter point for later (shape & syst lookup):
      self.param_of[(mass_str, actual_proc_name)] = params
      # add base
      self.base_of[actual_proc_name] = process.name
      # Register with CH and set shapes using cbCopy(mass_str, actual_proc_name, ...)
      add(params, mass_str, actual_proc_name)
    elif self.model.param_dependent_bkg:
      for signal_proc in self.processes.values():
        if signal_proc.is_signal:
          # Get model parameters for current parameter point
          params = signal_proc.params
          mass_str = str(int(params['mass'])) if self.is_resonant else "*"
          actual_proc_name = proc
          # record for uncertainties:
          self.param_of[(mass_str, actual_proc_name)] = params
          self.base_of[actual_proc_name] = proc
          add(params, mass_str, proc)
    else:
      # Backgrounds (and data) are stored under '*' unless bkg depends on params
      add(None, "*", proc)

  def addUncertainty(self, unc_name):
    unc = self.uncertainties[unc_name]
    isMVLnUnc = isinstance(unc, MultiValueLnNUncertainty)
    
    # Make sure param-independent backgrounds also get uncertainties
    items = list(self.param_of.items())
    if not self.model.param_dependent_bkg:
      for pname, p in self.processes.items():
        if p.is_background:
          items.append((( "*", pname ), None))  # params=None

    for (mass_str, process_name), params in items:
      for era, channel, category in self.ECC():
        #fall back for resonant
        base_name = self.base_of.get(process_name, process_name)
        process = self.processes[base_name]
        if process.is_data: continue
        
        # Get model parameters from param_of mapping
        model_params = params
        if isMVLnUnc:
          unc_value = self.getMultiValueLnUnc(unc,unc_name,process, era, channel, category, model_params)

        uncApplies = unc_value != None if isMVLnUnc else unc.appliesTo(process, era, channel, category)
        if not uncApplies: continue
        if not process.hasCompatibleModelParams(model_params, self.model.param_dependent_bkg): continue

        # Use mass_str and process_name from param_of mapping
        actual_proc_name = process_name
        actual_mass_str = mass_str

        nominal_shape = None
        shapes = {}
        if unc.needShapes:
          # Use model_params from param_of mapping
          nominal_shape = self.getShape(process, era, channel, category, model_params)
          for unc_scale in [ UncertaintyScale.Up, UncertaintyScale.Down ]:
            shapes[unc_scale] = self.getShape(process, era, channel, category, model_params,
                                              unc_name, unc_scale.name)
        unc_to_apply = unc.resolveType(nominal_shape, shapes, self.autolnNThr, self.asymlnNThr)
        can_ignore = unc_to_apply.canIgnore(unc_value, self.ignorelnNThr) if isMVLnUnc else unc_to_apply.canIgnore(self.ignorelnNThr)
        if can_ignore:
          print(f"Ignoring uncertainty {unc_name} for {process_name} in {era} {channel} {category}")
          continue
        systMap = unc_to_apply.valueToMap(unc_value) if isMVLnUnc else unc_to_apply.valueToMap()
        # Use actual process name and mass parameter
        cb_copy = self.cbCopy(actual_mass_str, actual_proc_name, era, channel, category)
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
          # Use actual process name and mass parameter
          cb_copy = self.cbCopy(actual_mass_str, actual_proc_name, era, channel, category).syst_name([unc_name])
          cb_copy.ForEachSyst(setShape)

  def writeDatacards(self, output):
    os.makedirs(output, exist_ok=True)
    
    # Simplified process collection - get all process names from CombineHarvester object (these now include unique signal names)
    all_process_names = set()
    for proc in self.cb.cp().process_set():
      all_process_names.add(proc)
    
    # Create separate datacards for each channel-category combination
    def slug(s: str) -> str:
      # lower, turn slashes into underscores, and collapse any weird chars to "_"
      s = s.lower().replace('/', '_')
      return re.sub(r'[^a-z0-9_.-]+', '_', s)
    
    for subera in self.eras:
      for subchannel in self.channels:
        for subcat in self.categories:
          channel_lower = subchannel.lower()               # eTau -> etau
          cat_slug = slug(subcat)                          # OS_Iso/res2b_cat3 -> os_iso_res2b_cat3
          dc_name = f"{channel_lower}_{cat_slug}"          # etau_os_iso_res2b_cat3

          # Output dir per (era, combo)
          tmp_output = os.path.join(output, subera, dc_name)
          os.makedirs(tmp_output, exist_ok=True)

          # Files
          tmp_dc_file = os.path.join(tmp_output, f"datacard_{dc_name}.txt")
          tmp_shape_file = os.path.join(tmp_output, f"shapes_{dc_name}.root")

          # Select exactly this bin (real CH name with slash is fine). mass was hardcoded as '*' for writing datacards
          bin_name = self.getBin(subera, subchannel, subcat, return_index=False)
          self.cb.cp().era([subera]).channel([subchannel]).bin([bin_name]).mass(['*']).WriteDatacard(
              tmp_dc_file, tmp_shape_file
          )

    # Write main combined datacard with all signals and backgrounds
    dc_file = os.path.join(output, "datacard_combined_signals.txt")
    # Create main combined shape file
    main_shape_file = os.path.join(output, "combined_signals_all.root")
    self.cb.cp().mass(['*']).WriteDatacard(dc_file, main_shape_file)



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
