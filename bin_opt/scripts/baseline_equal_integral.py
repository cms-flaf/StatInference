import ROOT, numpy, os, math, array

# the rebinAndFill function is from StatInference common tools:
def rebinAndFill(new_hist, old_hist, epsilon=1e-7):
  def check_range(old_axis, new_axis):
    old_min = old_axis.GetBinLowEdge(1)
    old_max = old_axis.GetBinUpEdge(old_axis.GetNbins())
    new_min = new_axis.GetBinLowEdge(1)
    new_max = new_axis.GetBinUpEdge(new_axis.GetNbins())
    return old_min <= new_min and old_max >= new_max

  def get_new_bin(old_axis, new_axis, bin_id_old):
    old_low_edge = round(old_axis.GetBinLowEdge(bin_id_old), 4)
    old_up_edge = round(old_axis.GetBinUpEdge(bin_id_old), 4)
    bin_low_new = new_axis.FindFixBin(old_low_edge)
    bin_up_new = new_axis.FindFixBin(old_up_edge)

    new_up_edge = new_axis.GetBinUpEdge(bin_low_new)
    if not (bin_low_new == bin_up_new or \
            abs(old_up_edge - new_up_edge) <= epsilon * abs(old_up_edge + new_up_edge) * 2):
      old_bins = [ str(old_axis.GetBinLowEdge(n)) for n in range(1, old_axis.GetNbins() + 2)]
      new_bins = [ str(new_axis.GetBinLowEdge(n)) for n in range(1, new_axis.GetNbins() + 2)]
      print('old_bins: [{}]'.format(', '.join(old_bins)))
      print('new_bins: [{}]'.format(', '.join(new_bins)))

      raise RuntimeError("Incompatible bin edges")
    return bin_low_new

  def add_bin_content(bin_old, bin_new):
    cnt_old = old_hist.GetBinContent(bin_old)
    err_old = old_hist.GetBinError(bin_old)
    cnt_new = new_hist.GetBinContent(bin_new)
    err_new = new_hist.GetBinError(bin_new)
    cnt_upd = cnt_new + cnt_old
    err_upd = math.sqrt(err_new ** 2 + err_old ** 2)
    new_hist.SetBinContent(bin_new, cnt_upd)
    new_hist.SetBinError(bin_new, err_upd);

  n_dim = old_hist.GetDimension()
  if new_hist.GetDimension() != n_dim:
    raise RuntimeError("Incompatible number of dimensions")
  if n_dim < 1 or n_dim > 2:
    raise RuntimeError("Unsupported number of dimensions")

  if not check_range(old_hist.GetXaxis(), new_hist.GetXaxis()):
    raise RuntimeError("x ranges are not compatible")

  if n_dim > 1 and not check_range(old_hist.GetYaxis(), new_hist.GetYaxis()):
    raise RuntimeError("y ranges are not compatible")

  for x_bin_old in range(1, old_hist.GetNbinsX() + 1):
    x_bin_new = get_new_bin(old_hist.GetXaxis(), new_hist.GetXaxis(), x_bin_old)
    if n_dim == 1:
      add_bin_content(x_bin_old, x_bin_new)
    else:
      for y_bin_old in range(1, old_hist.GetNbinsY() + 1):
        y_bin_new = get_new_bin(old_hist.GetYaxis(), new_hist.GetYaxis(), y_bin_old)
        bin_old = old_hist.GetBin(x_bin_old, y_bin_old)
        bin_new = new_hist.GetBin(x_bin_new, y_bin_new)
        add_bin_content(bin_old, bin_new)

def baseline_equal_integral(shape_file, num_bins, signal_processes):
    new_bin_edges = [0.0]
    new_bin_content = []
    new_bin_error = []

    for idx,signal in enumerate(signal_processes):
        if idx==0:
            signal_hist = shape_file.Get(signal)
            print(idx, signal, signal_hist.Integral())
        else:
            signal_hist.Add(shape_file.Get(signal))
            print(idx, signal, signal_hist.Integral())
    print("signal_hist integral()", signal_hist.Integral())

    total_signal_integral = signal_hist.Integral()
    print("total_signal_integral", total_signal_integral)

    target_integral_per_bin = total_signal_integral / num_bins
    print("target_integral_per_bin", target_integral_per_bin, "signal bins", signal_hist.GetNbinsX())

    tmp_integral = 0.0
    tmp_error2 = 0.0
    for bin in range(1, signal_hist.GetNbinsX()+1):
        bin_content = signal_hist.GetBinContent(bin)
        tmp_integral += bin_content
        tmp_error2 += signal_hist.GetBinError(bin)**2
        
        if tmp_integral >= target_integral_per_bin:
            edge = signal_hist.GetXaxis().GetBinUpEdge(bin)
            if edge > new_bin_edges[-1]:
                new_bin_edges.append(edge)
                new_bin_content.append(tmp_integral)
                new_bin_error.append(numpy.sqrt(tmp_error2))
                tmp_integral = 0.0
                tmp_error2 = 0.0
                print(bin, "edges:", new_bin_edges, "bin_content:", new_bin_content, "bin_error:", new_bin_error, "integral:", sum(new_bin_content))

    if new_bin_edges[-1] < 1.0:
        new_bin_edges.append(1.0)

    return new_bin_edges

def transform_binning(hist, binning, custom_bin_edges):
    def clone_hist(*args, **kwargs):
        clone = hist.__class__(
            hist.GetName()+'_clone',
            ";".join([hist.GetTitle(), hist.GetXaxis().GetTitle(), hist.GetYaxis().GetTitle()]),
            *args,
            **kwargs
        )
        clone.Sumw2()
        return clone

    if binning == "numbers":
        n_bins = hist.GetNbinsX()
        x_min = 0
        x_max = n_bins + x_min
        clone = clone_hist(n_bins, x_min, x_max)

        # fill values
        for b in range(1, n_bins + 1):
            clone.SetBinContent(b, hist.GetBinContent(b))
            clone.SetBinError(b, hist.GetBinError(b))
        for b in range(0, n_bins+1):
            clone.GetXaxis().ChangeLabel(b+1, labText=f"{custom_bin_edges[b]}")

        clone.GetXaxis().SetNdivisions(len(custom_bin_edges), 2, 2)

    elif binning == "numbers_width":
        n_bins = hist.GetNbinsX()
        x_min = 0.5
        x_max = n_bins + 0.5
        clone = clone_hist(n_bins, x_min, x_max)

        # fill values
        for b in range(1, n_bins + 1):
            clone.SetBinContent(b, hist.GetBinContent(b) / hist.GetBinWidth(b))
            clone.SetBinError(b, hist.GetBinError(b) / hist.GetBinWidth(b))

    return clone


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("--input-shape-file", required=True, type=str, help="input shape file")
    parser.add_argument("--num_bins", required=True, type=int, help="number of bins with equal integrals")
    parser.add_argument("--signal-process", required=True, type=str, help="list of signal process(es), i.e. ggHH_kl_1_kt_1_hbbhtt,ggHH_kl_2p45_kt_1_hbbhtt")

    args = parser.parse_args()

    shape_file = ROOT.TFile.Open(args.input_shape_file)
    num_bins = args.num_bins
    signal_process_list = [str(x) for x in args.signal_process.split(',')]
    # print(signal_process_list)

    output_dir = os.path.dirname(args.input_shape_file)
    # print(output_dir)
    if not os.path.exists(output_dir+"/equal_integral"):
        os.makedirs(output_dir+"/equal_integral")
    output_dir += "/equal_integral"
    input_file_name = os.path.basename(args.input_shape_file)
    # print(input_file_name)
    # print(os.path.splitext(input_file_name)[0])
    # print(os.path.join(output_dir, os.path.splitext(input_file_name)[0]+"_equal_integral_rebinned.root"))
    output_rebinned_file = ROOT.TFile(os.path.join(output_dir, os.path.splitext(input_file_name)[0]+"_equal_integral_rebinned_test.root"), "RECREATE")

    new_bin_edges = baseline_equal_integral(shape_file, num_bins, signal_process_list)
    print(new_bin_edges)

    for key in shape_file.GetListOfKeys():
        print("key.GetName()", key.GetName())
        hist = key.ReadObj()
        # print(hist.GetName(), hist.GetTitle())

        if hist.InheritsFrom("TH1"):
            rebinned = hist.Rebin(len(new_bin_edges)-1, hist.GetTitle(), array.array('d', new_bin_edges))
            # rebinned = ROOT.TH1D(hist.GetName(), hist.GetName(), len(new_bin_edges)-1, array.array('d', new_bin_edges))
            # rebinAndFill(rebinned, hist)
            rebinAndFill(rebinned, hist)
            rebinned_equidistant = transform_binning(rebinned, "numbers", custom_edges)

        else:
            continue

        # out_dir.cd()
        output_rebinned_file.cd()
        #rebinned.Write()
        rebinned_equidistant.Write()
        rebinned.Reset()
        rebinned_equidistant.Reset()

    output_rebinned_file.Close()
    shape_file.Close()

# python3 baseline_equal_integral.py --input-shape-file /eos/user/t/toakhter/bin_opt_tests/bin_opt_test_oct_hamburg_v2/hh_res2b_tauTau_2022_13p6TeV.input.root --num_bins 7 --signal-process ggHH_kl_1_kt_1_hbbhtt,ggHH_kl_2p45_kt_1_hbbhtt,ggHH_kl_5_kt_1_hbbhtt