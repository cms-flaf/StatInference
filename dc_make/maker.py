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
    CategoryNaming,
)
from .process import Process
from .uncertainty import (
    Uncertainty,
    UncertaintyType,
    UncertaintyScale,
    MultiValueLnNUncertainty,
    LnNUncertainty,
    ShapeUncertainty,
)
from .model import Model
from .binner import Binner

ROOT = importROOT()


class DatacardMaker:
    customizeble_parameters = ["eras", "channels", "categories"]

    def __init__(
        self,
        cfg_file,
        input_path,
        hist_bins=None,
        param_values=None,
        **kwargs,
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
        # The categories are exactly the ones the configuration lists -- the datacard bins
        # are whatever directories the input shapes contain, and nothing here derives them.
        # For input that a 2D->1D rebinning produced that means the sliced names
        # ("SR/res2b_dnn0"), and `category_pattern` is how they are taken apart again to
        # group the slices of one base category.
        self.naming = CategoryNaming.fromConfig(cfg)
        self.categories = list(cfg["categories"])
        self.signalFractionForRelevantBins = cfg["signalFractionForRelevantBins"]

        self.era_groups = cfg.get("era_groups", {})

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
        self._merged_away = {}
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

    def mergedAwayIn(self, channel, category):
        """Process names absorbed by an active merged process in this bin.

        A process declared with `subprocesses` that name other *datacard* processes
        replaces them wherever it applies -- see MinorBkg in the bbWW DL config,
        which merges DY/ST/VV in the boosted slices where DY alone has no usable MC
        statistics. Suppressing the constituents here is what stops the merge from
        double counting, and it means they keep their own configuration unchanged
        instead of needing a mirror-image category list to carve the merged bins
        back out.

        Inert for every existing configuration: the `subprocesses` lists in
        x_hh_bbtautau_run2.yaml name sample-level histograms (WW, WZ, ZZ ...), none
        of which is a datacard process, so nothing is ever absorbed there.
        """
        key = (channel, category)
        if key not in self._merged_away:
            absorbed = set()
            for p in self.processes.values():
                if not p.subprocesses:
                    continue
                if p.name not in self.channel_processes[channel]:
                    continue
                if not p.appliesToCategory(category):
                    continue
                absorbed |= {s for s in p.subprocesses if s in self.processes}
            self._merged_away[key] = absorbed
        return self._merged_away[key]

    def processInBin(self, name, channel, category):
        """Whether process `name` enters the datacard for this (channel, category)."""
        if name not in self.channel_processes[channel]:
            return False
        if not self.processes[name].appliesToCategory(category):
            return False
        return name not in self.mergedAwayIn(channel, category)

    def getSubEras(self, era):
        """Get sub-eras for a given era. If era is a meta-era, return its sub-eras.
        Otherwise return [era]."""
        if self.isMetaEra(era):
            return self.era_groups[era]
        return [era]

    def isMetaEra(self, era):
        """Check if era is a meta-era."""
        return era in self.era_groups

    def cbCopy(self, param_str, process, era, channel, category):
        bin_idx, bin_name = self.getBin(era, channel, category)
        return self.cb.cp().mass([param_str]).process([process]).bin([bin_name])

    def ECC(self):
        return itertools.product(self.eras, self.channels, self.categories)

    @staticmethod
    def lnNIsEraDependent(unc):
        """Whether an lnN's value depends on which era it is applied to.

        Only an era-dependent lnN needs the meta-era shape treatment. One that applies
        uniformly -- no `eras:` anywhere in its entry -- scales the whole summed yield by
        the same factor, which is what an lnN already says, so turning it into a template
        buys nothing and states the same uncertainty a different way: combine morphs a
        shape and takes an lnN as an exact log-normal.

        For MultiValueLnNUncertainty the eras are read off the value keys, not off
        `unc.eras`. Uncertainty.fromConfig reassigns its `args` dict inside the sub-entry
        loop and passes the result to the constructor, so a multi-value entry carries
        whatever the *last* sub-entry happened to scope by -- top_mass ends up with
        processes ('ST',) and lumi_1_13p6TeV with the 2023 eras. Nothing reads those today
        (addUncertainty decides multi-value applicability from getUncertaintyForProcess),
        but they are not to be trusted.
        """
        if isinstance(unc, MultiValueLnNUncertainty):
            # key is (processes, eras, channels, categories), built in fromConfig.
            return any(key[1] for key in unc.values)
        return bool(unc.eras)

    def uncAppliesTo(self, unc, process, era, channel, category):
        """Whether an uncertainty applies, with a meta-era inheriting its sub-eras'.

        Uncertainty.appliesTo matches the era it is given against the uncertainty's own
        `eras:` list, and a meta-era is not in it -- CMS_scale_j_2022 is scoped to
        Run3_2022, and Run3_Early is not that string. Answering False is right at that
        level; the meta-era is a fact this class knows and uncertainty.py does not.

        A nuisance scoped to one sub-era does apply to the combination: it varies that
        era's part of the summed shape. Without this every per-era shape systematic was
        dropped from a meta-era card silently -- 64 of 76 for HH->bbWW, including the jet
        energy scale and resolution. Per-era lnN escaped only because they are intercepted
        for meta-eras earlier, in _addMetaEraLnNAsShapeUnc.
        """
        eras = self.getSubEras(era) if self.isMetaEra(era) else [era]
        return any(unc.appliesTo(process, e, channel, category) for e in eras)

    def getCategoryGroups(self):
        """{"SR/res2b": ["SR/res2b_dnn0", ...]} -- the slices of one base category.

        Slices of the same base category are the natural unit for a per-category
        breakdown: they are one physical selection cut into pieces, not independent
        categories.
        """
        groups = {}
        for cat in self.categories:
            groups.setdefault(self.naming.base(cat), []).append(cat)
        return groups

    def PPECC(self):
        param_bins = list(self.param_bins.keys())
        if not self.model.param_dependent_bkg:
            param_bins.append("*")
        return itertools.product(
            self.processes.keys(), param_bins, self.eras, self.channels, self.categories
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

    @staticmethod
    def readHist(file, hist_name):
        """A histogram from a file, owned by Python rather than by ROOT.

        TFile.Get() hands back an object ROOT keeps alive for the lifetime of the file, so
        dropping the Python reference frees nothing. That is affordable when the objects
        are the small 1D shapes a datacard bin is made of; it is not when they are the 2D
        inputs a binning is cut from, where one era's build walked through ~15 GB. Taking
        ownership lets refcounting free each one after it has been sliced, while anything
        still referenced -- a cached shape -- stays alive as usual.

        Returns None if the path is missing: TFile.Get() on a fully-missing nested path
        returns a PyROOT wrapper around a null pointer, which is not `is None` but is
        falsy, so `if obj:` is the correct check.
        """
        obj = file.Get(hist_name)
        if not obj:
            return None
        ROOT.SetOwnership(obj, True)
        return obj

    def _loadBinnedHist(self, file, era, channel, category, model_params, hist_name):
        hist = self.readHist(file, hist_name)
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
        file_name, file = self.getInputFile(sub_era, model_params)
        hist_names = (
            [(subp, subp) for subp in process.subprocesses]
            if process.subprocesses
            else [(process.hist_name, process.name)]
        )
        up_hist = None
        down_hist = None
        applies = False

        for hist_name_suffix, proc_name_for_unc in hist_names:
            hist = self._loadBinnedHist(
                file,
                sub_era,
                channel,
                category,
                model_params,
                f"{channel}/{category}/{hist_name_suffix}",
            )
            unc_value = self._getLnNValue(
                unc, process, proc_name_for_unc, sub_era, channel, category
            )
            if unc_value is not None:
                applies = True
                sub_up = self._applyLnNToHist(hist, unc_value, UncertaintyScale.Up)
                sub_down = self._applyLnNToHist(hist, unc_value, UncertaintyScale.Down)
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

        sub_eras = self.getSubEras(era)

        if process.is_asimov_data:
            # Build the combined asimov sum from each background's own combined
            # (already negative-bin-resolved, with that background's own
            # tolerance) shape -- not from raw per-sub-era background shapes
            # summed then checked under data_obs's own (untolerant) settings.
            # This mirrors how a real era builds asimov data: by summing
            # already-resolved per-process shapes, never raw ones.
            combined_hist = None
            for bkg_proc in self.processes.values():
                if bkg_proc.is_background:
                    if not self.processInBin(bkg_proc.name, channel, category):
                        continue
                    bkg_hist = self.getCombinedShape(
                        bkg_proc, era, channel, category, model_params
                    )
                    if bkg_hist is None:
                        continue
                    if combined_hist is None:
                        combined_hist = bkg_hist.Clone()
                    else:
                        combined_hist.Add(bkg_hist)
            if combined_hist is None:
                raise RuntimeError("Cannot create asimov data histogram")
            return combined_hist

        # Meta-era: combine histograms from all sub-eras. Negative-bin
        # validation is deferred until after summing (below) rather than
        # applied per sub-era here -- a sub-era can dip negative on its own
        # statistics while the combined shape is fine, and only the combined
        # shape is what actually goes into the meta-era datacard.
        combined_hist = None

        for sub_era in sub_eras:
            # A per-era nuisance varies its own era and leaves the rest alone, so a sub-era
            # it is not scoped to contributes its *nominal*. Two things depend on that: the
            # variation only exists in its own era's file, so asking the others for it
            # raises; and the summed Up shape has to carry the full yield, or its ratio to
            # the nominal -- which is what the fit reads -- would be meaningless.
            applies = unc_name is None or self.uncAppliesTo(
                self.uncertainties[unc_name], process, sub_era, channel, category
            )
            sub_hist = self._rawShape(
                process,
                sub_era,
                channel,
                category,
                model_params,
                unc_name if applies else None,
                unc_scale if applies else None,
            )
            if sub_hist is None:
                continue
            if combined_hist is None:
                combined_hist = sub_hist.Clone()
            else:
                combined_hist.Add(sub_hist)

        needs_check = combined_hist is not None and not (
            process.is_signal and not (unc_name and unc_scale)
        )
        if needs_check:
            self.resolveOrRaiseNegativeBins(
                combined_hist,
                process,
                era,
                channel,
                category,
                model_params,
                unc_name,
                unc_scale,
                discovery_eras=sub_eras,
            )

        return combined_hist

    def _rawShape(
        self,
        process,
        era,
        channel,
        category,
        model_params,
        unc_name=None,
        unc_scale=None,
    ):
        """One process's shape in one bin, as read, with no negative-bin validation.

        Separate from getShape() for getCombinedShape()'s sake: a meta-era's sub-eras are
        summed before anything is checked, because a sub-era can dip negative on its own
        statistics while the combined shape -- the one that actually enters the datacard --
        is fine. Validating per sub-era would reject shapes that are good, and would also
        change them: resolveNegativeBins rebalances in place, so clamping each sub-era and
        then summing does not give the same shape as summing and then clamping.
        """
        file_name, file = self.getInputFile(era, model_params)
        key = (
            file_name,
            process.name,
            era,
            channel,
            category,
            unc_name,
            unc_scale,
        )

        if key not in self.shapes:
            if process.is_data and (unc_name is not None or unc_scale is not None):
                raise RuntimeError("Cannot apply uncertainty to the data process")
            if process.is_asimov_data:
                hist = None
                for bkg_proc in self.processes.values():
                    if bkg_proc.is_background:
                        if not self.processInBin(bkg_proc.name, channel, category):
                            continue
                        # Validated, not raw: asimov data is the sum of the
                        # per-process shapes as they enter the fit, each resolved
                        # under its own tolerance, never raw shapes summed and then
                        # checked under data_obs's own (untolerant) settings. The
                        # meta-era branch of getCombinedShape mirrors this.
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
                        subhist = self.readHist(file, hist_name)
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
                    hist = self.readHist(file, hist_name)
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
            self.shapes[key] = hist
        return self.shapes[key]

    def _isUnvalidatedSignal(self, process, unc_name, unc_scale):
        """Whether this is a nominal signal shape, which is never negative-bin checked.

        Those are collected into signal_hists_by_key to define the relevant (signal
        carrying) bins that the check itself consults, so they are the input to the rule
        rather than subject to it.
        """
        return process.is_signal and not (unc_name and unc_scale)

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
        """One process's shape in one bin, negative-bin validated -- what a datacard bin
        is built from.

        A meta-era is summed from its sub-eras first (getCombinedShape), which validates
        the sum rather than the parts; see _rawShape() for why.
        """
        if self.isMetaEra(era):
            return self.getCombinedShape(
                process, era, channel, category, model_params, unc_name, unc_scale
            )

        raw = self._rawShape(
            process, era, channel, category, model_params, unc_name, unc_scale
        )
        if raw is None or self._isUnvalidatedSignal(process, unc_name, unc_scale):
            return raw

        # Cached separately from the raw form, and validated on a copy: the check
        # rebalances in place, and the raw shape is still needed unmodified by
        # getCombinedShape, which sums the sub-eras before checking anything.
        key = (
            self.model.getInputFileName(era, model_params),
            process.name,
            era,
            channel,
            category,
            unc_name,
            unc_scale,
            "validated",
        )
        if key not in self.shapes:
            hist = raw.Clone()
            hist.SetDirectory(0)
            self.resolveOrRaiseNegativeBins(
                hist,
                process,
                era,
                channel,
                category,
                model_params,
                unc_name,
                unc_scale,
            )
            self.shapes[key] = hist
        return self.shapes[key]

    def resolveOrRaiseNegativeBins(
        self,
        hist,
        process,
        era,
        channel,
        category,
        model_params,
        unc_name=None,
        unc_scale=None,
        discovery_eras=None,
    ):
        """Validate/rebalance negative bins in-place on `hist`, raising if the
        result isn't accepted. `discovery_eras`, when given (meta-era combined
        shapes), unions relevant-signal-bin lookups across those real sub-eras
        instead of the single `era` -- signal shapes are cached per real
        sub-era, never under the meta-era name itself."""
        param_str = self.model.paramStr(model_params) if model_params else "*"
        key_param = (
            param_str
            if not self.keep_all_signal_hypothesis_into_single_datacard
            else "*"
        )
        lookup_eras = discovery_eras if discovery_eras else [era]
        signals = []
        for lookup_era in lookup_eras:
            signals.extend(
                self.signal_hists_by_key.get(
                    (lookup_era, channel, category, key_param), []
                )
            )
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

        final_integral = sum(
            hist.GetBinContent(n) for n in range(1, hist.GetNbinsX() + 1)
        )
        is_degenerate = not solution.accepted or final_integral <= 0

        if is_degenerate and unc_name and unc_scale:
            # A shape systematic variation that can't be resolved into a
            # valid (positive-integral) histogram -- whether flagged directly
            # by resolveNegativeBins, or only zero/negative after its donor
            # balancing happened to cancel out the whole shape -- is
            # inherently unusable for combine's shape
            # morphing (it requires a nonzero norm for every variation). This
            # is a low-statistics artifact of the up/down reweighting, not a
            # real central-value problem, so fall back to the nominal shape:
            # i.e. treat the systematic as having no effect in this bin.
            nominal = self.getShape(process, era, channel, category, model_params)
            for n in range(1, hist.GetNbinsX() + 1):
                hist.SetBinContent(n, nominal.GetBinContent(n))
                hist.SetBinError(n, nominal.GetBinError(n))
            return

        if not solution.accepted:
            axis = hist.GetXaxis()
            bins_edges = [
                str(axis.GetBinLowEdge(n)) for n in range(1, axis.GetNbins() + 2)
            ]
            bin_values = [
                str(hist.GetBinContent(n)) for n in range(1, axis.GetNbins() + 1)
            ]
            bin_errors = [
                str(hist.GetBinError(n)) for n in range(1, axis.GetNbins() + 1)
            ]
            print(f'bins_edges: [ {", ".join(bins_edges)} ]')
            print(f'bin_values: [ {", ".join(bin_values)} ]')
            print(f'bin_errors: [ {", ".join(bin_errors)} ]')
            raise RuntimeError(
                f"Negative bins found in histogram for {channel}/{category}/{process.hist_name}"
                + (f" (syst {unc_name}{unc_scale})" if unc_name and unc_scale else "")
            )

    def getSignalProcessForParams(self, model_params):
        """Signal Process matching model_params, or None. Used to gate a
        param-dependent background on whether the signal hypothesis it's
        being evaluated for actually has a shape in a given era/channel/
        category -- backgrounds are looked up per-MX (param_dependent_bkg),
        so a category the rebinning skipped for that MX has no background
        histograms either, not just no signal."""
        for p in self.processes.values():
            if p.is_signal and p.params == model_params:
                return p
        return None

    def hasNominalShape(self, process, era, channel, category):
        """Whether process's nominal shape exists for (era, channel, category),
        without raising. Used to skip a signal (and its per-mass background
        counterpart) where a specific era+category+mass genuinely has no signal
        MC -- e.g. a standalone single-era limit for a sparse category/channel
        that only has signal statistics once combined with other eras."""
        sub_eras = self.getSubEras(era) if self.isMetaEra(era) else [era]
        hist_name = f"{channel}/{category}/{process.hist_name}"
        for sub_era in sub_eras:
            _, file = self.getInputFile(sub_era, process.params)
            obj = self.readHist(file, hist_name)
            if obj is not None and obj.InheritsFrom("TH1"):
                return True
        return False

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
            if not self.hasNominalShape(process, era, channel, category):
                print(
                    f"Skipping {process.name} in {era}/{channel}/{category}: no signal shape found"
                )
                return
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
            # One copy of this process per distinct signal *parameter point*, not per
            # signal process: several signal processes (e.g. the bbWW and bbtautau
            # decay modes) share the same mass grid, and adding the copy once per
            # process would set the same shape twice ("Shape already set").
            seen_params = set()
            for signal_proc in self.processes.values():
                if not signal_proc.is_signal:
                    continue
                if not self.hasNominalShape(signal_proc, era, channel, category):
                    continue
                model_params = signal_proc.params
                param_str = (
                    self.model.paramStr(model_params)
                    if not self.keep_all_signal_hypothesis_into_single_datacard
                    else "*"
                )
                if param_str in seen_params:
                    continue
                seen_params.add(param_str)
                add(model_params, param_str, proc)
                self.param_of[(param_str, proc)] = model_params
                self.base_of[proc] = proc
        else:
            add(None, "*", proc)

    def addUncertainty(self, unc_name):
        unc = self.uncertainties[unc_name]
        isMVLnUnc = isinstance(unc, MultiValueLnNUncertainty)

        for proc, param_str, era, channel, category in self.PPECC():
            if not self.processInBin(proc, channel, category):
                continue
            process = self.processes[proc]
            if process.is_data:
                continue
            model_params = self.param_bins.get(param_str, None)
            if not process.hasCompatibleModelParams(
                model_params, self.model.param_dependent_bkg
            ):
                continue
            if process.is_signal:
                if not self.hasNominalShape(process, era, channel, category):
                    continue
            elif self.model.param_dependent_bkg and model_params is not None:
                signal_proc = self.getSignalProcessForParams(model_params)
                if signal_proc is not None and not self.hasNominalShape(
                    signal_proc, era, channel, category
                ):
                    continue

            if (
                self.isMetaEra(era)
                and isinstance(unc, (LnNUncertainty, MultiValueLnNUncertainty))
                and self.lnNIsEraDependent(unc)
            ):
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
                else self.uncAppliesTo(unc, process, era, channel, category)
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
            if isinstance(unc_to_apply, ShapeUncertainty):
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

                    if self.isMetaEra(era) and isinstance(
                        unc, (LnNUncertainty, MultiValueLnNUncertainty)
                    ):
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
                        else self.uncAppliesTo(unc, process, era, channel, category)
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

                    if isinstance(unc_to_apply, ShapeUncertainty):

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

        # Group the signal processes by parameter point. Several signal processes can
        # share a mass (e.g. the bbWW and bbtautau decay modes of the same resonance);
        # they must go into the *same* datacard so the fit scales them with a common
        # signal strength, rather than yielding a separate limit per decay mode.
        signals_by_param = {}
        for proc_name, process in self.processes.items():
            if not process.is_signal:
                continue
            key = self.model.paramStr(process.params)
            signals_by_param.setdefault(key, []).append(proc_name)

        for param_str, signal_names in signals_by_param.items():
            processes = list(signal_names) + background_names
            param_list = [param_str]
            if not self.model.param_dependent_bkg:
                param_list.append("*")
            # Named after the primary (first configured) signal, so a single-signal
            # config keeps exactly the file names it produced before.
            proc_name = signal_names[0]
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

                # Same breakdown by base category (all its slices, all channels),
                # for per-category limits alongside the per-channel ones.
                for base_cat, slice_cats in self.getCategoryGroups().items():
                    bin_names = [
                        self.getBin(subera, subchannel, cat, return_index=False)
                        for subchannel in self.channels
                        for cat in slice_cats
                    ]
                    selected = (
                        self.cb.cp()
                        .era([subera])
                        .bin(bin_names)
                        .mass(param_list)
                        .process(processes)
                    )
                    # A base category can be absent for a given mass hypothesis (e.g.
                    # boosted at low MX, where the rebinning found too little signal
                    # to slice it) -- there is no card to write then.
                    if len(selected.bin_set()) == 0:
                        continue
                    cat_dir = os.path.join(
                        output, subera, "categories", base_cat.replace("/", "_")
                    )
                    os.makedirs(cat_dir, exist_ok=True)
                    selected.WriteDatacard(
                        os.path.join(cat_dir, f"datacard_{proc_name}.txt"), shape_file
                    )

            self.cb.cp().mass(param_list).process(processes).WriteDatacard(
                dc_file, shape_file
            )

    def createDatacards(self, output, verbose=1):
        try:
            for era, channel, category in self.ECC():
                for name, p in self.processes.items():
                    if not self.processInBin(name, channel, category):
                        continue
                    if p.is_signal:
                        self.addProcess(name, era, channel, category)
            for era, channel, category in self.ECC():
                for name, p in self.processes.items():
                    if not self.processInBin(name, channel, category):
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
