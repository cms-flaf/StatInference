import json
from StatInference.common.tools import listToVector, rebinAndFill, importROOT
ROOT = importROOT()

class Binner:
    def __init__(self, hist_bins):
        self.hist_bins = []
        if type(hist_bins) == str:
            #Load json
            with open(hist_bins) as f:
                self.hist_bins = json.load(f)
        elif type(hist_bins) == list:
            self.hist_bins.append({
               'bins': hist_bins
            })           
        elif hist_bins is not None:
            raise RuntimeError("Incompatible hist_bins format")
        for entry in self.hist_bins:
           entry['bins'] = listToVector(entry['bins'], 'double')

    def applyBinning(self, era, channel, category, model_params, hist):
        if len(self.hist_bins) == 0:
            return hist
        new_binning = []
        def entry_passes(entry):
            if not ('eras' not in entry or era in entry['eras']): return False
            if not ('channels' not in entry or channel in entry['channels']): return False
            if not ('categories' not in entry or category in entry['categories']): return False

            for param_key, param_value in model_params.items():
                if not (param_key not in entry or param_value in entry[param_key]): return False
            return True

        for entry in self.hist_bins:
            if entry_passes(entry): new_binning.append(entry['bins'])
        if len(new_binning) <= 0:
            raise RuntimeError(f"No binning found for era/channel/category/params {era}/{channel}/{category}/{model_params}")
        if len(new_binning) >= 2:
            raise RuntimeError(f"Multiple binnings found for era/channel/category/params {era}/{channel}/{category}/{model_params}")

        new_hist = ROOT.TH1F(hist.GetName(), hist.GetTitle(), len(new_binning[0]) - 1,new_binning[0].data())
        rebinAndFill(new_hist, hist)
        return new_hist

