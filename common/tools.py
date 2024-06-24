import math

class PackageWrapper:
  def __init__(self, import_fn):
    self._package = None
    self._import_fn = import_fn

  def __getattr__(self, name):
    if self._package is None:
      self._package = self._import_fn()
    return getattr(self._package, name)

def lazyImport(package_name):
  def import_fn():
    return __import__(package_name)
  return PackageWrapper(import_fn)

def importROOT():
  def import_fn():
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.SetMustClean(False)
    ROOT.TH1.AddDirectory(False)
    ROOT.TH2.AddDirectory(False)
    ROOT.EnableThreadSafety()
    return ROOT
  return PackageWrapper(import_fn)

ROOT = importROOT()

def listToVector(x, value_type):
  v = ROOT.std.vector(value_type)(len(x))
  for n in range(len(x)):
    v[n] = x[n]
  return v

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

  for x_bin_old in range(old_hist.GetNbinsX() + 2):
    x_bin_new = get_new_bin(old_hist.GetXaxis(), new_hist.GetXaxis(), x_bin_old)
    if n_dim == 1:
      add_bin_content(x_bin_old, x_bin_new)
    else:
      for y_bin_old in range(old_hist.GetNbinsY() + 2):
        y_bin_new = get_new_bin(old_hist.GetYaxis(), new_hist.GetYaxis(), y_bin_old)
        bin_old = old_hist.GetBin(x_bin_old, y_bin_old)
        bin_new = new_hist.GetBin(x_bin_new, y_bin_new)
        add_bin_content(bin_old, bin_new)

def hasRelevantNegativeBins(histogram, relevant_bins = []):
  for n in range(1, histogram.GetNbinsX() + 1):
    if histogram.GetBinContent(n) < 0 and len(relevant_bins)>0 and relevant_bins[n] == True :
      print(f"bin {n} is relevant for signals and has negative content")
      return True
  return False

def hasNegativeBins(histogram):
  k=1
  for n in range(1, histogram.GetNbinsX() + 1):
    if histogram.GetBinContent(n) < 0 :
      if round(abs(histogram.GetBinContent(n)),10)+ k*round(abs(histogram.GetBinError(n)),10)>=0:
        return False
      print(f"for bin {n} histogram has {histogram.GetBinContent(n)} and error {histogram.GetBinError(n)}")
      return True
  return False
