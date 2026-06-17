import itertools
import math
import os
import yaml

from CombineHarvester.CombineTools.ch import CombineHarvester

from StatInference.common.tools import (
    listToVector,
    rebinAndFill,
    importROOT,
    resolveNegativeBins,
    getRelevantBins,
)
from .process import Process
from .uncertainty import (
    Uncertainty,
    UncertaintyType,
    UncertaintyScale,
    MultiValueLnNUncertainty,
    ShapeUncertainty,
)
from .model import Model
from .binner import Binner

ROOT = importROOT()


class DatacardMaker:
    customizeble_parameters = ["eras", "channels", "categories"]

    def __init__(
        self, cfg_file, input_path, hist_bins=None, param_values=None, **kwargs
    ):
        self.cb = CombineHarvester()

        self.input_path = input_path
        with open(cfg_file, "r") as f:
            cfg = yaml.safe_load(f)
        for param in self.customizeble_parameters:
            if param not in kwargs:
                continue
            value = kwargs[param]
            if value is None:
                continue
            value_type = type(cfg[param])
            if value_type == list:
                value = value.split(",")
            else:
                value = value_type(value)
            print(f"Overriding config parameter {param} to {value}")
            cfg[param] = value

        self.analysis = cfg["analysis"]
        self.eras = cfg["eras"]
        self.channels = cfg["channels"]
        self.categories = cfg["categories"]
        self.signalFractionForRelevantBins = cfg["signalFractionForRelevantBins"]

        # Load era groups if defined (e.g., Early_Run3 combining multiple eras)
        self.era_groups = cfg.get("era_groups", {})
        # Create reverse mapping: sub_era -> meta_era
        self.era_to_group = {}
        for group_name, sub_eras in self.era_groups.items():
            for sub_era in sub_eras:
                if sub_era not in self.era_to_group:
                    self.era_to_group[sub_era] = []
                self.era_to_group[sub_era].append(group_name)

        # Expand eras to include both regular eras and meta-eras
        all_eras = list(self.eras)
        for group_name in self.era_groups.keys():
            if group_name not in all_eras:
                all_eras.append(group_name)
        self.all_eras = all_eras

        self.bins = []
        for era, channel, cat in self.ECC():
            bin = self.getBin(era, channel, cat, return_index=False)
            self.bins.append(bin)

        self.model = Model.fromConfig(cfg["model"])
        self.keep_all_signal_hypothesis_into_single_datacard = cfg.get(
            "keep_all_signal_hypothesis_into_single_datacard", False
        )
        self.param_bins = {}
        self.processes = {}
        self.param_of = {}
        self.base_of = {}
        data_process = None
        has_signal = False
        self.channel_processes = {}
        for channel in self.channels:
            self.channel_processes[channel] = []

        for process in cfg["processes"]:
            if (type(process) != str) and process.get("is_signal", False):
                if param_values is not None:
                    print(f"Overwriting signal parameters to {param_values}")
                    process["param_values"] = param_values
            new_processes = Process.fromConfig(process, self.model)
            for process in new_processes:
                if process.name in self.processes:
                    raise RuntimeError(f"Process name {process.name} already exists")
                print(f"Adding {process}")
                self.processes[process.name] = process
                if process.channels:
                    for channel in process.channels:
                        if channel not in self.channel_processes:
                            print(f"Channel {channel} not defined in config")
                            continue
                        self.channel_processes[channel].append(process.name)
                else:
                    for channel in self.channels:
                        self.channel_processes[channel].append(process.name)
                if process.is_data:
                    if data_process is not None:
                        raise RuntimeError("Multiple data processes defined")
                    data_process = process
                if process.is_signal:
                    has_signal = True
                    self.param_bins.setdefault(
                        self.model.paramStr(process.params), process.params
                    )
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

        self.autoMCStats = cfg.get("autoMCStats", {"apply": False})

        hist_bins = hist_bins or cfg.get("hist_bins", None)
        self.hist_binner = Binner(hist_bins)
        # print(f"Using hist_bins: {self.hist_binner.hist_bins}")

        self.input_files = {}
        self.shapes = {}
        self.signal_hists_by_key = {}

    def getBin(self, era, channel, category, return_name=True, return_index=True):
        name = f"{era}_{self.analysis}_{channel}_{category}"
        if not return_name and not return_index:
            raise RuntimeError("Invalid argument combination")
        if not return_index:
            return name
        index = self.bins.index(name)
        if not return_name:
            return index
        return (index, name)

    def getSubEras(self, era):
        """Get sub-eras for a given era. If era is a meta-era, return its sub-eras.
        Otherwise return [era]."""
        if era in self.era_groups:
            return self.era_groups[era]
        return [era]

    def isMetaEra(self, era):
        """Check if era is a meta-era."""
        return era in self.era_groups

    def cbCopy(self, param_str, process, era, channel, category):
        bin_idx, bin_name = self.getBin(era, channel, category)
        return self.cb.cp().mass([param_str]).process([process]).bin([bin_name])

    def ECC(self):
        return itertools.product(self.all_eras, self.channels, self.categories)

    def PPECC(self):
        param_bins = list(self.param_bins.keys())
        if not self.model.param_dependent_bkg:
            param_bins.append("*")
        return itertools.product(
            self.processes.keys(),
            param_bins,
            self.all_eras,
            self.channels,
            self.categories,
        )

    def getInputFile(self, era, model_params):
        file_name = self.model.getInputFileName(era, model_params)
        if file_name not in self.input_files:
            full_file_name = os.path.join(self.input_path, file_name)
            file = ROOT.TFile.Open(full_file_name, "READ")
            if file == None:
                raise RuntimeError(f"Cannot open file {full_file_name}")
            self.input_files[file_name] = file
        return file_name, self.input_files[file_name]

    def _getLnNValue(self, unc, process, proc_name_for_unc, sub_era, channel, category):
        if isinstance(unc, MultiValueLnNUncertainty):
            return unc.getUncertaintyForProcess(
                proc_name_for_unc, sub_era, channel, category
            )
        if unc.appliesTo(process, sub_era, channel, category):
            return unc.value
        return None

    def _applyLnNToHist(self, hist, unc_value, direction):
        scaled = hist.Clone()
        if isinstance(unc_value, dict):
            factor = 1 + unc_value[direction]
        elif direction == UncertaintyScale.Up:
            factor = 1 + unc_value
        else:
            factor = 1 - unc_value
        scaled.Scale(factor)
        scaled.SetDirectory(0)
        return scaled

    def _loadBinnedHist(self, file, era, channel, category, model_params, hist_name):
        hist = file.Get(hist_name)
        if hist is None:
            raise RuntimeError(f"Cannot find histogram {hist_name} in {file.GetName()}")
        binned = self.hist_binner.applyBinning(
            era, channel, category, model_params, hist
        )
        binned.SetDirectory(0)
        return binned

    def _getSubEraLnNVariedShapes(
        self, unc, process, sub_era, channel, category, model_params
    ):
        if process.subprocesses and isinstance(unc, MultiValueLnNUncertainty):
            file_name, file = self.getInputFile(sub_era, model_params)
            up_hist = None
            down_hist = None
            applies = False

            for subp in process.subprocesses:
                hist = self._loadBinnedHist(
                    file,
                    sub_era,
                    channel,
                    category,
                    model_params,
                    f"{channel}/{category}/{subp}",
                )
                unc_value = self._getLnNValue(
                    unc, process, subp, sub_era, channel, category
                )
                if unc_value is not None:
                    applies = True
                    sub_up = self._applyLnNToHist(hist, unc_value, UncertaintyScale.Up)
                    sub_down = self._applyLnNToHist(
                        hist, unc_value, UncertaintyScale.Down
                    )
                else:
                    sub_up = hist.Clone()
                    sub_down = hist.Clone()
                    sub_up.SetDirectory(0)
                    sub_down.SetDirectory(0)

                if up_hist is None:
                    up_hist = sub_up
                    down_hist = sub_down
                else:
                    up_hist.Add(sub_up)
                    down_hist.Add(sub_down)

            if process.scale != 1:
                up_hist.Scale(process.scale)
                down_hist.Scale(process.scale)
            return up_hist, down_hist, applies

        nominal_sub = self.getShape(process, sub_era, channel, category, model_params)
        unc_value = self._getLnNValue(
            unc, process, process.name, sub_era, channel, category
        )
        if unc_value is None:
            up_hist = nominal_sub.Clone()
            down_hist = nominal_sub.Clone()
            up_hist.SetDirectory(0)
            down_hist.SetDirectory(0)
            return up_hist, down_hist, False

        up_hist = self._applyLnNToHist(nominal_sub, unc_value, UncertaintyScale.Up)
        down_hist = self._applyLnNToHist(nominal_sub, unc_value, UncertaintyScale.Down)
        return up_hist, down_hist, True

    def getMetaEraLnNShapeUnc(self, unc, process, era, channel, category, model_params):
        if not self.isMetaEra(era):
            return None

        nominal_shape = self.getShape(process, era, channel, category, model_params)
        combined_up = None
        combined_down = None
        any_applies = False

        for sub_era in self.getSubEras(era):
            up, down, applies = self._getSubEraLnNVariedShapes(
                unc, process, sub_era, channel, category, model_params
            )
            if applies:
                any_applies = True
            if combined_up is None:
                combined_up = up.Clone()
                combined_down = down.Clone()
            else:
                combined_up.Add(up)
                combined_down.Add(down)

        if not any_applies:
            return None
        return nominal_shape, {
            UncertaintyScale.Up: combined_up,
            UncertaintyScale.Down: combined_down,
        }

    def _canIgnoreLnNShape(self, nominal_shape, shapes):
        nom_int = nominal_shape.Integral()
        if nom_int == 0:
            return True
        up_frac = (shapes[UncertaintyScale.Up].Integral() - nom_int) / nom_int
        down_frac = (shapes[UncertaintyScale.Down].Integral() - nom_int) / nom_int
        return abs(up_frac) < self.ignorelnNThr and abs(down_frac) < self.ignorelnNThr

    def _addMetaEraLnNAsShapeUnc(
        self, unc_name, proc, param_str, process, era, channel, category, model_params
    ):
        unc = self.uncertainties[unc_name]
        shape_result = self.getMetaEraLnNShapeUnc(
            unc, process, era, channel, category, model_params
        )
        if shape_result is None:
            return False
        nominal_shape, shapes = shape_result
        if self._canIgnoreLnNShape(nominal_shape, shapes):
            print(
                f"Ignoring uncertainty {unc_name} for {proc} in {era} {channel} {category}"
            )
            return False

        cb_copy = self.cbCopy(param_str, proc, era, channel, category)
        cb_copy.AddSyst(
            self.cb,
            unc_name,
            UncertaintyType.shape.name,
            ShapeUncertainty(unc_name).valueToMap(),
        )
        shape_set = False

        def setShape(syst):
            nonlocal shape_set
            print(f"Setting unc shape for {syst}")
            if shape_set:
                raise RuntimeError("Shape already set")
            syst.set_shapes(
                shapes[UncertaintyScale.Up],
                shapes[UncertaintyScale.Down],
                nominal_shape,
            )
            shape_set = True

        self.cbCopy(param_str, proc, era, channel, category).syst_name(
            [unc_name]
        ).ForEachSyst(setShape)
        return True

    def getCombinedShape(
        self,
        process,
        era,
        channel,
        category,
        model_params,
        unc_name=None,
        unc_scale=None,
    ):
        """Combine histograms from multiple sub-eras for a meta-era.
        For meta-eras, this sums histograms from constituent sub-eras.
        For regular eras, delegates to getShape."""
        if not self.isMetaEra(era):
            # Regular era - just get the shape normally
            return self.getShape(
                process, era, channel, category, model_params, unc_name, unc_scale
            )

        # Meta-era: combine histograms from all sub-eras
        sub_eras = self.getSubEras(era)
        combined_hist = None

        for sub_era in sub_eras:
            sub_hist = self.getShape(
                process, sub_era, channel, category, model_params, unc_name, unc_scale
            )
            if sub_hist is None:
                continue
            if combined_hist is None:
                combined_hist = sub_hist.Clone()
            else:
                combined_hist.Add(sub_hist)

        return combined_hist

    def getShape(
        self,
        process,
        era,
        channel,
        category,
        model_params,
        unc_name=None,
        unc_scale=None,
    ):
        # Handle meta-eras by combining sub-era shapes
        if self.isMetaEra(era):
            return self.getCombinedShape(
                process, era, channel, category, model_params, unc_name, unc_scale
            )

        file_name, file = self.getInputFile(era, model_params)
        key = (file_name, process.name, era, channel, category, unc_name, unc_scale)

        if key not in self.shapes:
            if process.is_data and (unc_name is not None or unc_scale is not None):
                raise RuntimeError("Cannot apply uncertainty to the data process")
            if process.is_asimov_data:
                hist = None
                for bkg_proc in self.processes.values():
                    if bkg_proc.is_background:
                        if bkg_proc.name not in self.channel_processes[channel]:
                            continue
                        bkg_hist = self.getShape(
                            bkg_proc, era, channel, category, model_params
                        )
                        if hist is None:
                            hist = bkg_hist.Clone()
                        else:
                            hist.Add(bkg_hist)
                if hist is None:
                    raise RuntimeError("Cannot create asimov data histogram")
            else:
                hist_name = f"{channel}/{category}/{process.hist_name}"
                hists = []
                if process.subprocesses:
                    for subp in process.subprocesses:
                        hist_name = f"{channel}/{category}/{subp}"
                        if unc_name and unc_scale:
                            hist_name += f"_{unc_name}_{unc_scale}"
                        subhist = file.Get(hist_name)
                        if subhist == None:
                            raise RuntimeError(
                                f"Cannot find histogram {hist_name} in {file.GetName()}"
                            )
                        hists.append(
                            self.hist_binner.applyBinning(
                                era, channel, category, model_params, subhist
                            )
                        )
                else:
                    if unc_name and unc_scale:
                        hist_name += f"_{unc_name}_{unc_scale}"
                    hist = file.Get(hist_name)
                    if hist == None:
                        raise RuntimeError(
                            f"Cannot find histogram {hist_name} in {file.GetName()}"
                        )
                    hists.append(
                        self.hist_binner.applyBinning(
                            era, channel, category, model_params, hist
                        )
                    )
                if len(hists) == 0:
                    raise RuntimeError(f"hist list is empty for file {file.GetName()}")
                hist = hists[0]
                if len(hists) > 1:
                    for histy in hists[1:]:
                        hist.Add(histy)

                hist.SetDirectory(0)
                if process.scale != 1:
                    hist.Scale(process.scale)

                if process.is_signal and not (unc_name and unc_scale):
                    param_str = (
                        self.model.paramStr(model_params) if model_params else "*"
                    )
                    key_sig = (
                        era,
                        channel,
                        category,
                        (
                            param_str
                            if not self.keep_all_signal_hypothesis_into_single_datacard
                            else "*"
                        ),
                    )
                    self.signal_hists_by_key.setdefault(key_sig, []).append(hist)
                else:
                    param_str = (
                        self.model.paramStr(model_params) if model_params else "*"
                    )
                    key_sig = (
                        era,
                        channel,
                        category,
                        (
                            param_str
                            if not self.keep_all_signal_hypothesis_into_single_datacard
                            else "*"
                        ),
                    )
                    signals = self.signal_hists_by_key.get(key_sig, [])
                    relevant_bins = getRelevantBins(
                        era,
                        channel,
                        category,
                        signals,
                        self.signalFractionForRelevantBins,
                    )
                    solution = resolveNegativeBins(
                        hist,
                        relevant_bins=relevant_bins,
                        allow_zero_integral=process.allow_zero_integral,
                        allow_negative_bins_within_error=process.allow_negative_bins_within_error,
                        max_n_sigma_for_negative_bins=process.max_n_sigma_for_negative_bins,
                        allow_negative_integral=process.allow_negative_integral,
                    )

                    if not solution.accepted:
                        axis = hist.GetXaxis()
                        bins_edges = [
                            str(axis.GetBinLowEdge(n))
                            for n in range(1, axis.GetNbins() + 2)
                        ]
                        bin_values = [
                            str(hist.GetBinContent(n))
                            for n in range(1, axis.GetNbins() + 1)
                        ]
                        bin_errors = [
                            str(hist.GetBinError(n))
                            for n in range(1, axis.GetNbins() + 1)
                        ]
                        print(f'bins_edges: [ {", ".join(bins_edges)} ]')
                        print(f'bin_values: [ {", ".join(bin_values)} ]')
                        print(f'bin_errors: [ {", ".join(bin_errors)} ]')
                        raise RuntimeError(
                            f"Negative bins found in histogram for {channel}/{category}/{process.hist_name}"
                            + (
                                f" (syst {unc_name}{unc_scale})"
                                if unc_name and unc_scale
                                else ""
                            )
                        )
            self.shapes[key] = hist
        return self.shapes[key]

    def addProcess(self, proc, era, channel, category):
        bin_idx, bin_name = self.getBin(era, channel, category)
        process = self.processes[proc]

        def add(model_params, param_str, process_name):
            if process.is_data:
                self.cb.AddObservations(
                    [param_str],
                    [self.analysis],
                    [era],
                    [channel],
                    [(bin_idx, bin_name)],
                )
            else:
                self.cb.AddProcesses(
                    [param_str],
                    [self.analysis],
                    [era],
                    [channel],
                    [process_name],
                    [(bin_idx, bin_name)],
                    process.is_signal,
                )

            shape = self.getShape(process, era, channel, category, model_params)
            shape_set = False

            def setShape(p):
                nonlocal shape_set
                print(f"Setting shape for {p}")
                if shape_set:
                    raise RuntimeError("Shape already set")
                p.set_shape(shape, True)
                shape_set = True

            cb_copy = self.cbCopy(param_str, process_name, era, channel, category)
            if process.is_data:
                cb_copy.ForEachObs(setShape)
            else:
                cb_copy.ForEachProc(setShape)

        if process.is_signal:
            model_params = process.params
            param_str = self.model.paramStr(model_params)
            if self.keep_all_signal_hypothesis_into_single_datacard:
                actual_proc_name = f"{process.name}_{param_str}"
                add(model_params, "*", actual_proc_name)
                self.param_of[("*", actual_proc_name)] = model_params
                self.base_of[actual_proc_name] = process.name
            else:
                actual_proc_name = process.name
                add(model_params, param_str, actual_proc_name)
                self.param_of[(param_str, actual_proc_name)] = model_params
                self.base_of[actual_proc_name] = process.name

        elif self.model.param_dependent_bkg:
            for signal_proc in self.processes.values():
                if not signal_proc.is_signal:
                    continue
                model_params = signal_proc.params
                param_str = (
                    self.model.paramStr(model_params)
                    if not self.keep_all_signal_hypothesis_into_single_datacard
                    else "*"
                )
                add(model_params, param_str, proc)
                self.param_of[(param_str, proc)] = model_params
                self.base_of[proc] = proc
        else:
            add(None, "*", proc)

    def addUncertainty(self, unc_name):
        unc = self.uncertainties[unc_name]
        isMVLnUnc = isinstance(unc, MultiValueLnNUncertainty)

        for proc, param_str, era, channel, category in self.PPECC():
            if proc not in self.channel_processes[channel]:
                continue
            process = self.processes[proc]
            if process.is_data:
                continue
            model_params = self.param_bins.get(param_str, None)
            if not process.hasCompatibleModelParams(
                model_params, self.model.param_dependent_bkg
            ):
                continue

            if self.isMetaEra(era) and unc.type == UncertaintyType.lnN:
                self._addMetaEraLnNAsShapeUnc(
                    unc_name,
                    proc,
                    param_str,
                    process,
                    era,
                    channel,
                    category,
                    model_params,
                )
                continue

            if isMVLnUnc:
                unc_value = unc.getUncertaintyForProcess(
                    process.name, era, channel, category
                )

            uncApplies = (
                unc_value != None
                if isMVLnUnc
                else unc.appliesTo(process, era, channel, category)
            )
            if not uncApplies:
                continue

            nominal_shape = None
            shapes = {}
            if unc.needShapes:
                model_params = self.param_bins.get(param_str, None)
                nominal_shape = self.getShape(
                    self.processes[proc], era, channel, category, model_params
                )

                for unc_scale in [UncertaintyScale.Up, UncertaintyScale.Down]:
                    shapes[unc_scale] = self.getShape(
                        self.processes[proc],
                        era,
                        channel,
                        category,
                        model_params,
                        unc_name,
                        unc_scale.name,
                    )

            unc_to_apply = unc.resolveType(
                nominal_shape, shapes, self.autolnNThr, self.asymlnNThr
            )
            can_ignore = (
                unc_to_apply.canIgnore(unc_value, self.ignorelnNThr)
                if isMVLnUnc
                else unc_to_apply.canIgnore(self.ignorelnNThr)
            )
            if can_ignore:
                print(
                    f"Ignoring uncertainty {unc_name} for {proc} in {era} {channel} {category}"
                )
                continue
            systMap = (
                unc_to_apply.valueToMap(unc_value)
                if isMVLnUnc
                else unc_to_apply.valueToMap()
            )
            cb_copy = self.cbCopy(param_str, proc, era, channel, category)
            cb_copy.AddSyst(self.cb, unc_name, unc_to_apply.type.name, systMap)
            if unc_to_apply.type == UncertaintyType.shape:
                shape_set = False

                def setShape(syst):
                    nonlocal shape_set
                    print(f"Setting unc shape for {syst}")
                    if shape_set:
                        raise RuntimeError("Shape already set")
                    syst.set_shapes(
                        shapes[UncertaintyScale.Up],
                        shapes[UncertaintyScale.Down],
                        nominal_shape,
                    )
                    shape_set = True

                self.cbCopy(param_str, proc, era, channel, category).syst_name(
                    [unc_name]
                ).ForEachSyst(setShape)

        if self.keep_all_signal_hypothesis_into_single_datacard:
            for (param_str, proc_name), params in self.param_of.items():
                for era, channel, category in self.ECC():
                    base_name = self.base_of.get(proc_name, proc_name)
                    process = self.processes[base_name]
                    if process.is_data:
                        continue
                    if not process.hasCompatibleModelParams(
                        params, self.model.param_dependent_bkg
                    ):
                        continue

                    if self.isMetaEra(era) and unc.type == UncertaintyType.lnN:
                        self._addMetaEraLnNAsShapeUnc(
                            unc_name,
                            proc_name,
                            param_str,
                            process,
                            era,
                            channel,
                            category,
                            params,
                        )
                        continue

                    if isMVLnUnc:
                        unc_value = unc.getUncertaintyForProcess(
                            process.name, era, channel, category
                        )
                    uncApplies = (
                        (unc_value is not None)
                        if isMVLnUnc
                        else unc.appliesTo(process, era, channel, category)
                    )
                    if not uncApplies:
                        continue

                    nominal_shape = None
                    shapes = {}
                    if unc.needShapes:
                        nominal_shape = self.getShape(
                            process, era, channel, category, params
                        )
                        for us in [UncertaintyScale.Up, UncertaintyScale.Down]:
                            shapes[us] = self.getShape(
                                process,
                                era,
                                channel,
                                category,
                                params,
                                unc_name,
                                us.name,
                            )

                    unc_to_apply = unc.resolveType(
                        nominal_shape, shapes, self.autolnNThr, self.asymlnNThr
                    )
                    if (
                        isMVLnUnc
                        and unc_to_apply.canIgnore(unc_value, self.ignorelnNThr)
                    ) or (not isMVLnUnc and unc_to_apply.canIgnore(self.ignorelnNThr)):
                        continue

                    systMap = (
                        unc_to_apply.valueToMap(unc_value)
                        if isMVLnUnc
                        else unc_to_apply.valueToMap()
                    )
                    cb_copy = self.cbCopy(param_str, proc_name, era, channel, category)
                    cb_copy.AddSyst(self.cb, unc_name, unc_to_apply.type.name, systMap)

                    if unc_to_apply.type == UncertaintyType.shape:

                        def setShape(syst):
                            syst.set_shapes(
                                shapes[UncertaintyScale.Up],
                                shapes[UncertaintyScale.Down],
                                nominal_shape,
                            )

                        self.cbCopy(
                            param_str, proc_name, era, channel, category
                        ).syst_name([unc_name]).ForEachSyst(setShape)

    def writeDatacards(self, output):
        os.makedirs(output, exist_ok=True)

        if self.keep_all_signal_hypothesis_into_single_datacard:
            dc_file = os.path.join(output, "datacard_combined_signals.txt")
            main_shape_file = os.path.join(output, "combined_signals_all.root")
            self.cb.cp().mass(["*"]).WriteDatacard(dc_file, main_shape_file)

            for subera, subchannel, subcat in self.ECC():
                tmp_dir = os.path.join(output, subera, subchannel, subcat)
                os.makedirs(tmp_dir, exist_ok=True)
                _, bin_name = self.getBin(subera, subchannel, subcat)
                perbin_dc = os.path.join(tmp_dir, f"datacard_{bin_name}.txt")
                perbin_root = os.path.join(tmp_dir, f"{bin_name}.root")
                self.cb.cp().bin([bin_name]).mass(["*"]).WriteDatacard(
                    perbin_dc, perbin_root
                )

            return

        background_names = [n for n, p in self.processes.items() if p.is_background]
        for proc_name, process in self.processes.items():
            if not process.is_signal:
                continue
            processes = [proc_name] + background_names
            param_list = [self.model.paramStr(process.params)]
            if not self.model.param_dependent_bkg:
                param_list.append("*")
            dc_file = os.path.join(output, f"datacard_{proc_name}.txt")
            shape_file = os.path.join(output, f"{proc_name}.root")

            for subera in self.eras:
                for subchannel in self.channels:
                    tmp_output = os.path.join(output, subera, subchannel)
                    os.makedirs(tmp_output, exist_ok=True)
                    tmp_dc_file = os.path.join(tmp_output, f"datacard_{proc_name}.txt")
                    tmp_shape_file = shape_file
                    self.cb.cp().era([subera]).channel([subchannel]).mass(
                        param_list
                    ).process(processes).WriteDatacard(tmp_dc_file, tmp_shape_file)

            # Also write datacards for meta-eras
            for meta_era in self.era_groups.keys():
                for subchannel in self.channels:
                    tmp_output = os.path.join(output, meta_era, subchannel)
                    os.makedirs(tmp_output, exist_ok=True)
                    tmp_dc_file = os.path.join(tmp_output, f"datacard_{proc_name}.txt")
                    tmp_shape_file = shape_file
                    self.cb.cp().era([meta_era]).channel([subchannel]).mass(
                        param_list
                    ).process(processes).WriteDatacard(tmp_dc_file, tmp_shape_file)

            self.cb.cp().mass(param_list).process(processes).WriteDatacard(
                dc_file, shape_file
            )

    def createDatacards(self, output, verbose=1):
        try:
            for era, channel, category in self.ECC():
                for name, p in self.processes.items():
                    if name not in self.channel_processes[channel]:
                        continue
                    if p.is_signal:
                        self.addProcess(name, era, channel, category)
            for era, channel, category in self.ECC():
                for name, p in self.processes.items():
                    if name not in self.channel_processes[channel]:
                        continue
                    if not p.is_signal:
                        self.addProcess(name, era, channel, category)

            for unc_name in self.uncertainties.keys():
                print(f"adding uncertainty: {unc_name}")
                self.addUncertainty(unc_name)
            if self.autoMCStats["apply"]:
                self.cb.SetAutoMCStats(
                    self.cb,
                    self.autoMCStats["threshold"],
                    self.autoMCStats["apply_to_signal"],
                    self.autoMCStats["mode"],
                )
            if verbose > 0:
                self.cb.PrintAll()
            self.writeDatacards(output)
        finally:
            for file in self.input_files.values():
                file.Close()
