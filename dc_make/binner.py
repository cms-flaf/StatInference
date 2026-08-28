import array
import json
from StatInference.common.tools import (
    listToVector,
    rebinAndFill,
    importROOT,
    CategoryNaming,
)
import FLAF.Common.HistHelper as HistHelper

ROOT = importROOT()


class Binner:
    def __init__(self, hist_bins):
        self.hist_bins = []
        self.linearize = False
        # A binning derived by bin_opt_2d/rebin_2d.py: it cuts a 2D input into per-slice
        # 1D shapes, which the other two formats cannot express. Recognised by its own
        # keys rather than by a flag, so a configuration just points hist_bins at it.
        self.slices = None
        self.naming = None
        if type(hist_bins) == str:
            # Load json
            with open(hist_bins) as f:
                self.hist_bins = json.load(f)
            if (
                isinstance(self.hist_bins, dict)
                and "category_pattern" in self.hist_bins
            ):
                self.slices = self.hist_bins
                self.naming = CategoryNaming(self.slices["category_pattern"])
                self.hist_bins = []
        elif type(hist_bins) == list:
            self.hist_bins.append({"bins": hist_bins})
        elif hist_bins is not None:
            raise RuntimeError("Incompatible hist_bins format")

        # Loaded the string, now what is the object?
        for entry in self.hist_bins:
            if "combined_bins" in entry:
                self.linearize = True
            else:
                entry["bins"] = listToVector(entry["bins"], "double")

    def sliceOf(self, channel, category, model_params):
        """The recorded slice for one datacard bin, or None if there is no slice binning.

        The category is a sliced name ("SR/res2b_dnn0"); the record is keyed by the base
        category it was cut from, with the slices in index order. Both come from the one
        category_pattern the record carries, so the names cannot be read back differently
        than they were written.
        """
        if self.slices is None:
            return None
        base, slice_idx = self.naming.split(category)
        if slice_idx is None:
            raise RuntimeError(
                f"{category} is not a sliced category name for pattern "
                f"'{self.naming.pattern}', but the binning describes slices"
            )
        param = str(list(model_params.values())[0]) if model_params else None
        entry = self.slices["binning"].get(param, {}).get(channel, {}).get(base)
        if entry is None:
            return None
        return entry["slices"][slice_idx]

    def applySlice(self, hist2d, entry):
        """One slice of a 2D input as the 1D shape that becomes a datacard bin.

        x is integrated over the slice's own range and y is rebinned onto the recorded
        edges; the result is the same histogram the binning was derived to produce.
        """
        edges = array.array("d", entry["y_edges"])
        xlo, xhi = entry["x_range"]
        out = ROOT.TH1D(hist2d.GetName(), hist2d.GetTitle(), len(edges) - 1, edges)
        out.SetDirectory(0)
        for bin_idx, (ylo, yhi) in enumerate(entry["y_ranges"], start=1):
            err = array.array("d", [0.0])
            out.SetBinContent(bin_idx, hist2d.IntegralAndError(xlo, xhi, ylo, yhi, err))
            out.SetBinError(bin_idx, err[0])
        return out

    def applyBinning(self, era, channel, category, model_params, hist):
        if self.slices is not None:
            entry = self.sliceOf(channel, category, model_params)
            if entry is None:
                raise RuntimeError(
                    f"No binning recorded for {channel}/{category} at "
                    f"{model_params}; the datacard configuration asks for a bin the "
                    "rebinning did not produce"
                )
            return self.applySlice(hist, entry)
        if len(self.hist_bins) == 0:
            return hist
        new_binning = []

        def entry_passes(entry):
            if not ("eras" not in entry or era in entry["eras"]):
                return False
            if not ("channels" not in entry or channel in entry["channels"]):
                return False
            if not ("categories" not in entry or category in entry["categories"]):
                return False
            if model_params is not None:
                for param_key, param_value in model_params.items():
                    if not (param_key not in entry or param_value in entry[param_key]):
                        return False
            return True

        if self.linearize:
            for entry in self.hist_bins:
                if entry_passes(entry):
                    new_binning.append(entry)
            new_hist = HistHelper.RebinHisto(hist, new_binning[0], hist.GetTitle())
        else:
            for entry in self.hist_bins:
                # print(entry)
                if entry_passes(entry):
                    new_binning.append(entry["bins"])
            if len(new_binning) <= 0:
                raise RuntimeError(
                    f"No binning found for era/channel/category/params {era}/{channel}/{category}/{model_params}"
                )
            if len(new_binning) >= 2:
                raise RuntimeError(
                    f"Multiple binnings found for era/channel/category/params {era}/{channel}/{category}/{model_params}"
                )
            new_hist = ROOT.TH1F(
                hist.GetName(),
                hist.GetTitle(),
                len(new_binning[0]) - 1,
                new_binning[0].data(),
            )
            rebinAndFill(new_hist, hist)
        return new_hist
