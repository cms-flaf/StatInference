import argparse
import ROOT
import os
import math

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)


def reading_shape_file_histograms(input_hists):
    histogram_dicts = {}

    for key in input_hists.GetListOfKeys():
        if key.GetName().startswith("data"):
            if "data" not in histogram_dicts:
                histogram_dicts["data"] = {}
            histogram_dicts["data"][key.GetName()] = input_hists.Get(key.GetName())
            print('data', key.GetName(), input_hists.Get(key.GetName()).Integral())
        elif (key.GetName().startswith("ggHH") or key.GetName().startswith("qqHH")):
            if "signal" not in histogram_dicts:
                histogram_dicts["signal"] = {}
            histogram_dicts["signal"][key.GetName()] = input_hists.Get(key.GetName())
        elif not (key.GetName().startswith("data") or key.GetName().startswith("ggHH") or key.GetName().startswith("qqHH")):
            if "background" not in histogram_dicts:
                histogram_dicts["background"] = {}
            histogram_dicts["background"][key.GetName()] = input_hists.Get(key.GetName())
        else:
            raise RuntimeError("No background/signal/data histograms found with the matching strings. Update the directory in main and/or histogram keys strings in the function reading_shape_file_histograms")

    # file.Close()
    return histogram_dicts

def signal_hists(signal_hist_dict, sig_names):
    signal_dict = {}
    for sig in sig_names:
        signal_dict[sig] = {}
        signal_dict[sig]["hists"] = []
        signal_dict[sig]["colors"] = []
    ggF = None
    VBF = None
    for hist_name in signal_hist_dict:
        hist = signal_hist_dict[hist_name]
        if hist.GetName().startswith("ggHH"):
            if ggF is None:
                ggF = hist.Clone("ggF")
                ggF.SetDirectory(0)
                if not ggF.GetSumw2():
                    ggF.Sumw2()
            else:
                ggF.Add(hist)
                if not ggF.GetSumw2():
                    ggF.Sumw2()


        elif hist.GetName().startswith("qqHH"):
            if VBF is None:
                VBF = hist.Clone("VBF")
                VBF.SetDirectory(0)
                if not VBF.GetSumw2():
                    VBF.Sumw2()
            else:
                VBF.Add(hist)
                if not VBF.GetSumw2():
                    VBF.Sumw2()
            

    signal_dict["ggF"]["hists"].append(ggF)
    signal_dict["VBF"]["hists"].append(VBF)
    if signal_dict["ggF"]["colors"]==[]: signal_dict["ggF"]["colors"].append(ROOT.kMagenta-10)
    if signal_dict["VBF"]["colors"]==[]: signal_dict["VBF"]["colors"].append(ROOT.kTeal)
    
    print(signal_dict)
    return signal_dict


def background_hists(background_hist_dict, bkg_names):
    background_dict = {}
    for bkg in bkg_names:
        background_dict[bkg] = {}
        background_dict[bkg]["hists"] = []
        background_dict[bkg]["colors"] = []
    background_dict["Others"] = {}
    background_dict["Others"]["hists"] = []
    background_dict["Others"]["colors"] = []
    

    singleH = None
    others = None

    for hist_name in background_hist_dict:
        hist = background_hist_dict[hist_name]
        if hist.GetName()=="DY":
            background_dict["DY"]["hists"].append(hist)
            if background_dict["DY"]["colors"]==[]: background_dict["DY"]["colors"].append(ROOT.kSpring-5)

        elif hist.GetName()=="TT":
            background_dict["TT"]["hists"].append(hist)
            for n in range(1, hist.GetNbinsX()+1):
                print("TT bin:", n, "content:", hist.GetBinContent(n), "error:", hist.GetBinError(n))
            if background_dict["TT"]["colors"]==[]: background_dict["TT"]["colors"].append(ROOT.kOrange-2)

        elif hist.GetName()=="QCD":
            background_dict["QCD"]["hists"].append(hist)
            if background_dict["QCD"]["colors"]==[]: background_dict["QCD"]["colors"].append(ROOT.kOrange+1)

        elif hist.GetName() in ["WH_htt", "ZH_hbb", "ggH_htt", "qqH_htt", "ttH_hbb"]:
            if singleH is None:
                singleH = hist.Clone("singleH")
                singleH.SetDirectory(0)
                if not singleH.GetSumw2():
                    singleH.Sumw2()
            else:
                singleH.Add(hist)
                if not singleH.GetSumw2():
                    singleH.Sumw2()
            if background_dict["singleH"]["colors"]==[]: background_dict["singleH"]["colors"].append(ROOT.kAzure+10)

        else:
            if others is None:
                others = hist.Clone("others")
                others.SetDirectory(0)
                if not others.GetSumw2():
                    others.Sumw2()
            else:
                others.Add(hist)
                if not others.GetSumw2():
                    others.Sumw2()
            if background_dict["Others"]["colors"]==[]: background_dict["Others"]["colors"].append(ROOT.kRed-4)

    background_dict["singleH"]["hists"].append(singleH)
    background_dict["Others"]["hists"].append(others)

    total_bkg = None
    for bkg in background_dict:
        for hist in background_dict[bkg]["hists"]:

            if total_bkg is None:
                total_bkg = hist.Clone("total_bkg")
                total_bkg.SetDirectory(0)
                if not total_bkg.GetSumw2():
                    total_bkg.Sumw2()
            else:
                total_bkg.Add(hist)
                if not total_bkg.GetSumw2():
                    total_bkg.Sumw2()

    return background_dict, total_bkg

def uncertainty(total_bkg_hist):
    band = total_bkg_hist.Clone("bkg_rel_unc")
    band.SetDirectory(0)
    n_bins = band.GetNbinsX()
    if not band.GetSumw2():
        band.Sumw2()

    for i in range(1, n_bins + 1):
        bin_content = total_bkg_hist.GetBinContent(i)
        bin_error = total_bkg_hist.GetBinError(i)
        if bin_content > 0:
            band.SetBinContent(i, bin_content)
            band.SetBinError(i, bin_error)
    
    band.SetFillColor(ROOT.kGray+2)
    band.SetFillStyle(3004)
    band.SetLineColor(ROOT.kGray+2)

    return band

def plot_stack_ratio(data_hist, bkg_hists_dict, total_bkg_hist, signal_hist, uncertainty_band_hist, output_path, want_data):
    canvas = ROOT.TCanvas("canvas", "", 800, 800)

    pad_main = ROOT.TPad("pad_main", "Main", 0, 0.3, 1, 1)
    pad_ratio = ROOT.TPad("pad_ratio", "Ratio", 0, 0, 1, 0.3)
    pad_main.SetBottomMargin(0)
    pad_ratio.SetTopMargin(0)
    pad_ratio.SetBottomMargin(0.3)
    pad_main.SetLogy()
    pad_main.Draw()
    pad_ratio.Draw()

    pad_main.cd()
    stack = ROOT.THStack("stack", "")

    sorting_bkg_by_name = ["Others", "singleH", "QCD", "TT", "DY"] #["DY", "TT", "QCD", "singleH", "Others"]
    sorted_bkg_hists_dict = sorted(bkg_hist_dict.items(), key=lambda x: sorting_bkg_by_name.index(x[0]))
    for bkg in sorted_bkg_hists_dict:
        for hist in bkg[1]["hists"]:
            hist.SetFillColor(bkg[1]["colors"][0])
            # hist.SetLineColor(ROOT.kBlack)
            stack.Add(hist)

    stack.Draw("HIST")
    stack.GetYaxis().SetTitleOffset(1.2)
    stack.GetYaxis().SetTitle("Events")
    stack.SetMinimum(0.001) #0.01*stack.GetMinimum())
    stack.SetMaximum(1000)

    uncertainty_band_hist.Draw("E2 SAME")

    for sig in signal_hist:
        for hist in signal_hist[sig]["hists"]:
            hist.SetLineColor(signal_hist[sig]["colors"][0])
            hist.SetLineWidth(3)
            hist.Draw("HIST SAME")

    pad_ratio.cd()

    # if want_data:
    #     if data_hist:
    #         data_hist.SetMarkerStyle(20)
    #         data_hist.SetMarkerSize(1)
    #         data_hist.SetLineColor(ROOT.kBlack)
    #         #last bin blinded
    #         data_hist.SetBinContent(data_hist.GetNbinsX(), 0)
    #         data_hist.SetBinError(data_hist.GetNbinsX(), 0)

    #         data_hist.Draw("E SAME")


    #     ratio = data_hist.Clone("ratio")
    #     ratio.SetDirectory(0)
    #     if not ratio.GetSumw2():
    #         ratio.Sumw2()
    #     if not total_bkg_hist.GetSumw2():
    #         total_bkg_hist.Sumw2()
    #     ratio.Divide(total_bkg_hist)
    #     # for i in range(1, ratio.GetNbinsX()+1):
    #     #     data_content = data_hist.GetBinContent(i)
    #     #     bkg_content = total_bkg_hist.GetBinContent(i)
    #     #     data_error = data_hist.GetBinError(i)
    #     #     bkg_error = total_bkg_hist.GetBinError(i)

    #     #     if bkg_content != 0:
    #     #         ratio_content = data_content / bkg_content
    #     #         ratio_error = math.sqrt( (data_error / bkg_content)**2 + (data_content * bkg_error / (bkg_content**2))**2 )
    #     #     else:
    #     #         ratio_content = 0
    #     #         ratio_error = 0

    #     #     #last bin blinded
    #     #     if i==ratio.GetNbinsX():
    #     #         ratio_content = 0
    #     #         ratio_error = 2

    #     #     ratio.SetBinContent(i, ratio_content)
    #     #     ratio.SetBinError(i, ratio_error)
    #     #last bin blinded
    #     ratio.SetBinContent(ratio.GetNbinsX(), 0)
    #     ratio.SetBinError(ratio.GetNbinsX(), 2)

    #     ratio.GetYaxis().SetNdivisions(505)
    #     ratio.GetYaxis().SetRangeUser(0.0, 2.0)
    #     ratio.GetYaxis().SetTitle("Data / Bkg")
    #     ratio.GetYaxis().SetTitleSize(0.12)
    #     ratio.GetYaxis().SetLabelSize(0.08)

    #     ratio.GetXaxis().SetTitle("DNN HH output node")
    #     ratio.GetXaxis().SetTitleSize(0.12)
    #     ratio.GetXaxis().SetLabelSize(0.1)
        
    #     ratio.Draw("E")

    uncertainty_band_ratio = uncertainty_band_hist.Clone("bkg_rel_unc_ratio")
    for bin in range(1, uncertainty_band_ratio.GetNbinsX()+1):
        content = uncertainty_band_ratio.GetBinContent(bin)
        error = uncertainty_band_ratio.GetBinError(bin)
        if content != 0:
            rel_error = error / content
        else:
            rel_error = 0
        uncertainty_band_ratio.SetBinContent(bin, 1.0)
        uncertainty_band_ratio.SetBinError(bin, rel_error)

    uncertainty_band_ratio.SetDirectory(0)

    uncertainty_band_ratio.GetYaxis().SetNdivisions(505)
    uncertainty_band_ratio.GetYaxis().SetRangeUser(0.0, 2.0)
    uncertainty_band_ratio.GetYaxis().SetTitle("Data / Bkg")
    uncertainty_band_ratio.GetYaxis().SetTitleSize(0.12)
    uncertainty_band_ratio.GetYaxis().SetLabelSize(0.08)

    uncertainty_band_ratio.GetXaxis().SetTitle("DNN HH output node")
    uncertainty_band_ratio.GetXaxis().SetTitleSize(0.12)
    uncertainty_band_ratio.GetXaxis().SetLabelSize(0.1)

    if want_data:
        uncertainty_band_ratio.Draw("E2 SAME")
    else:
        uncertainty_band_ratio.Draw("E2")

    pad_main.cd()

    legend = ROOT.TLegend(0.65, 0.6, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    if data_hist:
        legend.AddEntry(data_hist, "Data", "ep")
    for bkg in bkg_hists_dict:
        for hist in bkg_hists_dict[bkg]["hists"]:
            legend.AddEntry(hist, hist.GetName(), "f")
    for sig in signal_hist:
        for hist in signal_hist[sig]["hists"]:
            legend.AddEntry(hist, "Expected "+sig+" HH", "l")
    legend.AddEntry(uncertainty_band_hist, "Bkg. uncert.", "f")
    
    legend.Draw()

    text_plot = ROOT.TLatex()
    text_plot.SetNDC()
    text_plot.SetTextSize(0.04)

    text_plot.SetTextFont(61)
    text_plot.DrawLatex(0.12, 0.92, "CMS")
    text_plot.SetTextFont(52)
    text_plot.DrawLatex(0.19, 0.92, "Preliminary")
    text_plot.SetTextFont(42)
    text_plot.DrawLatex(0.65, 0.92, "7.9804 fb^{-1} (13.6 TeV)")

    text_plot.SetTextSize(0.04)
    text_plot.DrawLatex(0.12, 0.87, "HH #rightarrow bb#tau#tau")
    text_plot.DrawLatex(0.25, 0.87, "res2b, #tau_{h}#tau_{h}")
    
    canvas.Update()
    canvas.SaveAs(output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True, type=str, help="input file directory and file name")
    parser.add_argument("--output", required=True, type=str, help="output directory and file name (.png/.pdf). if no directory given, will save in the cwd")
    parser.add_argument("--wantData", required=True, type=bool, help="whether to plot data or not")
    parser.add_argument("--signal-name", required=False, type=str, default="ggF,VBF", help="comma separated signal names")
    parser.add_argument("--bkg-names", required=True, type=str, help="comma-sparated names of background")

    args = parser.parse_args()

    input_file = ROOT.TFile.Open(args.input_file, "READ")
    file_dir = input_file.Get("tauTau/OS_Iso/res2b")

    hist_dict = reading_shape_file_histograms(file_dir)

    bkg_hist_dict, total_bkg_hist = background_hists(hist_dict["background"], args.bkg_names.split(","))
    signal_hist_dict = signal_hists(hist_dict["signal"], args.signal_name.split(","))
    uncertainty_band_hist = uncertainty(total_bkg_hist)
    plot_stack_ratio(hist_dict["data"]["data_obs"], bkg_hist_dict, total_bkg_hist, signal_hist_dict, uncertainty_band_hist, args.output, args.wantData)

    input_file.Close()

# $ python3 stackPlot.py --input_file /eos/user/t/toakhter/bin_opt_tests/bin_opt_test_oct_hamburg_v2/plotting_test/hh_res2b_tauTau_2022_13p6TeV.input.root --output /eos/user/t/toakhter/bin_opt_tests/bin_opt_test_oct_hamburg_v2/plotting_test/test_1.pdf --bkg-names DY,TT,QCD,singleH --wantData False
