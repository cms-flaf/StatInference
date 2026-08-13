import array
import math
import os
import sys
import yaml

if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    pkg_dir = os.path.dirname(file_dir)
    base_dir = os.path.dirname(pkg_dir)
    pkg_dir_name = os.path.split(pkg_dir)[1]
    if base_dir not in sys.path:
        sys.path.append(base_dir)
    __package__ = pkg_dir_name

from StatInference.common.tools import importROOT, CategoryNaming
from StatInference.common.param_parse import extractParameters, applyParameters
from StatInference.dc_make.model import Model

ROOT = importROOT()


def load_config(config_path):
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)
    model = Model.fromConfig(cfg["model"])
    channels = cfg["channels"]
    # cfg["categories"] holds base category names ("SR/res2b"), which is exactly what
    # the raw 2D histograms are keyed by. The per-slice names this script writes are
    # built by the CategoryNaming run() installs below.
    categories = list(dict.fromkeys(cfg["categories"]))

    # Several processes may carry is_signal (e.g. the bbWW(2l) and bbtautau decay
    # modes of the same resonance, both scaled by the same signal strength). They
    # are summed to form the discovery signal that steers the DNN slice
    # boundaries, so the binning is optimised for the total signal in the fit.
    signal_hist_names = []
    mass_values = None
    background_entries = []
    for entry in cfg["processes"]:
        if type(entry) == str:
            background_entries.append((entry, entry, []))
            continue
        if entry.get("is_data", False):
            continue
        base_name = entry["process"]
        hist_name = entry.get("hist_name", base_name)
        if entry.get("is_signal", False):
            signal_hist_names.append(hist_name)
            if mass_values is None:
                mass_values = entry["param_values"]
            elif list(entry["param_values"]) != list(mass_values):
                raise RuntimeError(
                    f"Signal {hist_name} has param_values {entry['param_values']}, "
                    f"which differ from {mass_values}; every signal must be defined "
                    "at the same mass points"
                )
        else:
            background_entries.append((base_name, hist_name, entry.get("channels", [])))

    if not signal_hist_names:
        raise RuntimeError("No signal process found in config")

    return {
        "model": model,
        "channels": channels,
        "categories": categories,
        "signal_hist_name_patterns": signal_hist_names,
        "signal_param_name": extractParameters(signal_hist_names[0])[0],
        "mass_values": mass_values,
        "background_entries": background_entries,
    }


def open_input_file(input_dir, model, era, mass, param_name):
    file_name = model.getInputFileName(era, {param_name: mass})
    full_path = os.path.join(input_dir, file_name)
    f = ROOT.TFile.Open(full_path, "READ")
    if f is None or f.IsZombie():
        raise RuntimeError(f"Cannot open file {full_path}")
    return f


def _detach(h):
    """Detach a histogram from its TFile and hand ownership to Python.

    SetDirectory(0) alone makes the histogram survive the file's close, but
    leaves it owned by nobody -- so every 2D histogram read here (one per
    systematic variation, per category, per mass) leaked for the lifetime of the
    process. Across the ten-mass loop that grew without bound and got the job
    SIGKILLed partway through the final mass, leaving a truncated ROOT file.
    SetOwnership makes the object die with its last Python reference.
    """
    h.SetDirectory(0)
    ROOT.SetOwnership(h, True)
    return h


def get_hist(f, path):
    h = f.Get(path)
    # TFile.Get() on a fully-missing nested path can return a PyROOT wrapper
    # around a null C++ pointer, which is not `is None` but is falsy.
    if not h:
        return None
    return _detach(h)


def sum_hists(hists):
    total = None
    for h in hists:
        if h is None:
            continue
        if total is None:
            total = _detach(h.Clone())
        else:
            total.Add(h)
    return total


def _integral(hist, lo, hi):
    return (
        hist.Integral(lo, hi, 0, -1)
        if hist.GetDimension() == 2
        else hist.Integral(lo, hi)
    )


def _integral_and_error(hist, lo, hi):
    """(yield, MC statistical error) over [lo, hi]."""
    err = array.array("d", [0.0])
    if hist.GetDimension() == 2:
        value = hist.IntegralAndError(lo, hi, 0, -1, err)
    else:
        value = hist.IntegralAndError(lo, hi, err)
    return value, err[0]


def _bkg_yields(bkg_hists_by_name, lo, hi):
    """For each background, the summed yield across all discovery eras (combined
    statistics). Individual eras are not required to individually clear the
    threshold -- use allow_negative_bins_within_error (maker.py, per process/
    category) for backgrounds/categories where a specific era can land negative
    in a bin that's fine once combined."""
    return {
        name: sum(_integral(h, lo, hi) for h in hists)
        for name, hists in bkg_hists_by_name.items()
    }


def _bkg_errors(bkg_hists_by_name, lo, hi):
    """Per background, the MC statistical error on the summed yield (eras added in
    quadrature)."""
    return {
        name: math.sqrt(sum(_integral_and_error(h, lo, hi)[1] ** 2 for h in hists))
        for name, hists in bkg_hists_by_name.items()
    }


def _total_bkg_error(bkg_hists_by_name, lo, hi):
    errors = _bkg_errors(bkg_hists_by_name, lo, hi)
    return math.sqrt(sum(e**2 for e in errors.values()))


def effective_entries(value, error):
    """(yield / error)^2 -- the unweighted MC event count a weighted yield is worth.

    A DY slice holding a couple of very-high-weight aMC@NLO events can carry a
    sizeable yield with an effective count far below 1, i.e. a background estimate
    that is statistically compatible with almost anything.
    """
    if error <= 0:
        return float("inf") if value > 0 else 0.0
    return (value / error) ** 2


def asimov_significance(s, b, b_err=0.0):
    """Median discovery significance for a counting experiment (Cowan et al.,
    arXiv:1007.1727), with the background uncertainty folded in when given.

    Reduces to sqrt(2*((s+b)*ln(1+s/b) - s)) for b_err = 0. Unlike S/sqrt(B) this
    stays valid when b is O(1), which is exactly the regime the high-mass boosted
    slices live in -- there S/sqrt(B) reports significances of ~30 on ~1.5
    background events, and drives the slice boundary on that basis.
    """
    if b <= 0 or s <= 0:
        return 0.0
    if b_err <= 0:
        return math.sqrt(max(2.0 * ((s + b) * math.log1p(s / b) - s), 0.0))
    var = b_err**2
    term1 = (s + b) * math.log(((s + b) * (b + var)) / (b * b + (s + b) * var))
    term2 = (b * b / var) * math.log1p(var * s / (b * (b + var)))
    return math.sqrt(max(2.0 * (term1 - term2), 0.0))


def significance(s, b, b_err=0.0, mode="sb"):
    """Figure of merit steering the DNN slice boundaries.

    mode="sb":     S / sqrt(B + sigma_B^2). Folding sigma_B into plain S/sqrt(B)
                   stops a downward fluctuation of a statistics-starved background
                   (DY in the b-tagged muMu slices, effective MC count <1) from
                   inflating the apparent significance. Still assumes the Gaussian
                   regime, which breaks down for B of order a few.
    mode="asimov": the Poisson-correct Asimov significance, valid at low B.
    """
    if b is None or b <= 0:
        return 0.0
    if mode == "asimov":
        return asimov_significance(s, b, b_err)
    denominator = b + b_err**2
    if denominator <= 0:
        return 0.0
    return s / (denominator**0.5)


def _dnn_passes(
    yields,
    min_sum,
    total_error=None,
    min_neff=0.0,
    errors=None,
    min_each=0.0,
    min_proc_neff=0.0,
    exempt=(),
):
    """DNN-slice validity: the summed background must clear min_sum, and be known
    to better than min_neff effective MC entries -- a window whose background is
    statistically undetermined must not be selectable at all.

    The per-process arms (min_each, min_proc_neff) are what _mass_passes already
    does on the mass axis, applied here as well. Testing only the sum hides an
    individual background that has fluctuated negative behind a large, well
    measured neighbour: muMu/SR/res2b at m500 selected a slice holding
    DY = -5.87 +- 8.92 (N_eff 0.43) because the summed background there was
    +39.9 +- 9.0 (N_eff 22), clearing both summed tests comfortably. That slice
    cannot be turned into a datacard at all -- a negative DY integral is rejected
    outright by resolveNegativeBins -- so the binning has to not select it in the
    first place. Worse, the significance being maximised is S/sqrt(B + sigma_B^2),
    which a downward fluctuation *increases* by lowering B, so such a window is
    mildly preferred rather than merely tolerated.

    `exempt` comes from minor_backgrounds() judged over the whole DNN range of the
    category, never over the candidate window: a background that is negligible in
    this category must not be able to veto every boundary, and judging it inside
    the window under test is circular (a background that is exactly zero there is
    trivially below any fraction of the total).
    """
    total = sum(yields.values())
    if total <= min_sum:
        return False
    if min_neff > 0 and total_error is not None:
        if effective_entries(total, total_error) < min_neff:
            return False
    if min_each > 0 or min_proc_neff > 0:
        for name, value in yields.items():
            if name in exempt:
                continue
            if value <= min_each:
                return False
            if min_proc_neff > 0 and errors is not None:
                if effective_entries(value, errors.get(name, 0.0)) < min_proc_neff:
                    return False
    return True


def minor_backgrounds(bkg_hists_by_name, lo, hi, min_frac):
    """Backgrounds negligible over the *whole* [lo, hi] range, which may therefore
    be exempted from the per-bin min_each floor.

    This must be judged once over the full slice, never inside the candidate bin
    under test. Evaluating the fraction per bin is circular: a background that is
    exactly zero in that bin is trivially below any fraction of the bin total, so
    it is always exempted -- which is precisely the case the floor exists to
    catch. That let a bin through with TT = 0 in a slice where TT is 94% of the
    background (m300, eMu, res2b_dnn2).
    """
    if min_frac <= 0:
        return set()
    yields = _bkg_yields(bkg_hists_by_name, lo, hi)
    total = sum(yields.values())
    if total <= 0:
        return set()
    return {name for name, value in yields.items() if value < min_frac * total}


def _mass_passes(yields, min_each, exempt=(), total_error=None, min_neff=0.0):
    """Mass-bin validity: every *relevant* background must exceed min_each, and
    (when min_neff > 0) the summed background must be known to at least min_neff
    effective MC entries.

    `exempt` is the set of backgrounds judged negligible across the whole slice by
    minor_backgrounds(); they are excused from min_each so that a process worth
    0.4% of the yield -- and known only to a few hundred percent -- cannot veto
    every candidate split, which would make find_mass_bins() back off all the way
    to a single bin and discard the mass shape in exactly the high-significance
    slices that matter most.

    The min_neff arm is the same gate _dnn_passes() applies to slice boundaries.
    Applying it only there left the mass bins *inside* each slice ungated, and
    since that is where essentially every fit bin lives, the median per-bin
    effective background count came out at ~2.8 against a slice threshold of 4.
    Note it constrains only the *summed* background, so it is satisfied by any one
    well-measured process; min_each/exempt is what protects the individual ones.
    """
    total = sum(yields.values())
    if min_neff > 0 and total_error is not None:
        if effective_entries(total, total_error) < min_neff:
            return False
    for name, value in yields.items():
        if name in exempt:
            continue
        if value <= min_each:
            return False
    return True


def grow_dnn_slice(
    sig_hist,
    bkg_hists_by_name,
    right,
    first_bin,
    min_sum,
    min_neff=0.0,
    sig_mode="sb",
    min_each=0.0,
    min_proc_neff=0.0,
    exempt=(),
):
    """Among every candidate [left, right] with summed background > min_sum,
    pick the one maximizing S/sqrt(B + sigma_B^2) -- not just the first one that
    clears it. This is what actually drives where the cut lands; the content
    floor is only a validity gate. Falls back to first_bin (best effort) if no
    candidate clears the floor anywhere."""
    best_left = None
    best_sig = -1.0
    for left in range(right, first_bin - 1, -1):
        bkg_y = _bkg_yields(bkg_hists_by_name, left, right)
        b_err = _total_bkg_error(bkg_hists_by_name, left, right)
        bkg_e = (
            _bkg_errors(bkg_hists_by_name, left, right) if min_proc_neff > 0 else None
        )
        if not _dnn_passes(
            bkg_y,
            min_sum,
            b_err,
            min_neff,
            bkg_e,
            min_each,
            min_proc_neff,
            exempt,
        ):
            continue
        s = _integral(sig_hist, left, right)
        b = sum(bkg_y.values())
        sig = significance(s, b, b_err, sig_mode)
        if sig > best_sig:
            best_sig = sig
            best_left = left
    return best_left if best_left is not None else first_bin


def find_dnn_slices(
    sig_hist,
    bkg_hists_by_name,
    n_slices,
    first_bin,
    last_bin,
    min_sum,
    min_neff=0.0,
    sig_mode="sb",
    min_each=0.0,
    min_proc_neff=0.0,
    min_frac=0.0,
):
    """Split [first_bin, last_bin] into exactly `n_slices` ranges (fixed count,
    required so every mass point shares the same category list), scanning right
    to left. Each slice (except the final, leftover one) is placed to maximize
    signal significance among boundaries with summed background > min_sum.

    Which backgrounds are minor enough to be exempt from the per-process floors is
    decided once here, over the whole [first_bin, last_bin] range, and held fixed
    for every candidate window -- see _dnn_passes on why it cannot be re-judged
    per window.
    """
    exempt = (
        minor_backgrounds(bkg_hists_by_name, first_bin, last_bin, min_frac)
        if (min_each > 0 or min_proc_neff > 0)
        else set()
    )
    slices = []
    right = last_bin
    for slice_idx in range(n_slices):
        if slice_idx == n_slices - 1 or right <= first_bin:
            lo = first_bin if right >= first_bin else right
            slices.append((lo, right if right >= lo else lo))
            right = lo - 1
            continue
        left = grow_dnn_slice(
            sig_hist,
            bkg_hists_by_name,
            right,
            first_bin,
            min_sum,
            min_neff,
            sig_mode,
            min_each,
            min_proc_neff,
            exempt,
        )
        slices.append((left, right))
        right = left - 1
    slices.reverse()
    return slices


def signal_quantile_ranges(sig_hist, n_bins, first_bin, last_bin):
    """Split [first_bin, last_bin] into `n_bins` ranges each holding an equal share
    of the signal.

    This is what puts the bins where the resonance is. The HME axis is 0-1500 GeV
    for every mass point, so the signal occupies a narrow window whose position
    moves with MX while the axis does not; binning must follow the signal rather
    than the axis. Equal-signal quantiles do that automatically -- bin edges
    cluster wherever dS/dm is large (the peak) and a single wide bin absorbs the
    long empty stretches on either side.

    The CDF clamps negative bin contents to zero so it stays monotonic; a
    statistical undershoot in a signal MC bin must not move an edge backwards.
    """
    n_avail = last_bin - first_bin + 1
    n_bins = max(1, min(n_bins, n_avail))
    if n_bins == 1:
        return [(first_bin, last_bin)]

    cumulative = []
    running = 0.0
    for b in range(first_bin, last_bin + 1):
        running += max(sig_hist.GetBinContent(b), 0.0)
        cumulative.append(running)
    total = running
    if total <= 0:
        return [(first_bin, last_bin)]

    ranges = []
    lo = first_bin
    for k in range(1, n_bins):
        target = k * total / n_bins
        b = lo
        while b < last_bin and cumulative[b - first_bin] < target:
            b += 1
        # every one of the n_bins-k ranges still to come needs at least one bin
        b = max(lo, min(b, last_bin - (n_bins - k)))
        ranges.append((lo, b))
        lo = b + 1
    ranges.append((lo, last_bin))
    return ranges


def merge_until_valid(
    ranges, sig_hist, bkg_hists_by_name, min_each, exempt=(), min_neff=0.0
):
    """Merge adjacent ranges until every one satisfies the background gates.

    The signal quantiles decide where the edges want to be; this decides how many
    of them the background statistics can actually support. A failing range is
    merged into whichever neighbour holds *less* signal, so the dense bins around
    the peak -- the ones carrying the discrimination -- are the last to be given
    up. Terminates because every step removes one range, ending at the single
    full-range bin, which has nothing left to fail against.
    """
    ranges = list(ranges)
    while len(ranges) > 1:
        bad = None
        for i, (lo, hi) in enumerate(ranges):
            if not _mass_passes(
                _bkg_yields(bkg_hists_by_name, lo, hi),
                min_each,
                exempt,
                _total_bkg_error(bkg_hists_by_name, lo, hi),
                min_neff,
            ):
                bad = i
                break
        if bad is None:
            break
        if bad == 0:
            other = 1
        elif bad == len(ranges) - 1:
            other = bad - 1
        else:
            left_sig = sig_hist.Integral(*ranges[bad - 1])
            right_sig = sig_hist.Integral(*ranges[bad + 1])
            other = bad - 1 if left_sig <= right_sig else bad + 1
        first, second = min(bad, other), max(bad, other)
        ranges[first : second + 1] = [(ranges[first][0], ranges[second][1])]
    return ranges


def find_mass_bins(
    sig_hist,
    bkg_hists_by_name,
    max_bins,
    first_bin,
    last_bin,
    min_each,
    min_frac=0.0,
    min_neff=0.0,
):
    """Mass bins inside one DNN slice: signal quantiles for the edges, background
    gates for the count.

    The previous version scanned right-to-left growing each bin leftward until the
    backgrounds cleared their floors, mirroring find_dnn_slices(). That is right
    for the DNN axis (signal-like = high score, so resolution belongs at the top)
    and wrong for HME, where the signal sits at HME ~ MX and both tails are empty.
    Starting from the top of the axis spent the bin budget on the empty
    high-HME tail and left the entire resonance peak in the single leftover bin:
    at MX=600 in muMu/res2b_dnn3, one bin covered 0-710 GeV holding 97.5% of the
    signal while five bins shared the 710-1500 GeV region holding 2.5%. The fit
    then had no shape to work with in the only region where signal and background
    differ.
    """
    # which backgrounds count as negligible is decided once, over the whole slice
    exempt = minor_backgrounds(bkg_hists_by_name, first_bin, last_bin, min_frac)
    ranges = signal_quantile_ranges(sig_hist, max_bins, first_bin, last_bin)
    return merge_until_valid(
        ranges, sig_hist, bkg_hists_by_name, min_each, exempt, min_neff
    )


def extend_outer_edges(ranges, full_lo, full_hi):
    """Widen the first/last range to swallow under/overflow (bin 0 / nbins+1),
    so no events are silently dropped at the extremes of the axis."""
    ranges = list(ranges)
    ranges[0] = (full_lo, ranges[0][1])
    ranges[-1] = (ranges[-1][0], full_hi)
    return ranges


def mass_bin_budget(bkg_hists_by_name, lo, hi, max_mass_bins, bkg_per_mass_bin):
    """How many mass bins this DNN slice can actually afford.

    A fixed max_mass_bins is applied blind to slice content: the high-mass boosted
    slices hold ~1.5 total background events and were still being split into 10
    bins, i.e. ~0.15 events per bin. Capping at B_slice / bkg_per_mass_bin ties the
    binning to the statistics that are really there. Returns max_mass_bins
    unchanged when bkg_per_mass_bin <= 0 (feature off).
    """
    if bkg_per_mass_bin <= 0:
        return max_mass_bins
    total = sum(_bkg_yields(bkg_hists_by_name, lo, hi).values())
    if total <= 0:
        return 1
    return max(1, min(max_mass_bins, int(total / bkg_per_mass_bin)))


def discover_binning(
    sig2d,
    bkg2d_by_name,
    n_dnn_slices,
    max_mass_bins,
    min_dnn_sum,
    min_mass_each,
    min_bkg_neff=0.0,
    min_bkg_frac=0.0,
    min_mass_bkg_neff=0.0,
    bkg_per_mass_bin=0.0,
    sig_mode="sb",
    min_dnn_bkg_each=0.0,
    min_dnn_bkg_neff=0.0,
):
    """bkg2d_by_name: {background_name: [hist per discovery era, ...]}. The list
    is usually a single era's own histogram (standalone limit) or all of a
    meta-era's sub-eras (combined limit) -- see HistRebinTask.get_discovery_eras().
    Yields are summed across whatever's in the list; see _bkg_yields(). sig2d is
    the same discovery reference's (already-summed) signal histogram: its x
    projection picks the significance-maximizing DNN slice boundaries, and its y
    projection within each slice places the mass bin edges by signal quantile."""
    any_hist = next(iter(bkg2d_by_name.values()))[0]
    nx = any_hist.GetNbinsX()
    ny = any_hist.GetNbinsY()
    dnn_slices = find_dnn_slices(
        sig2d,
        bkg2d_by_name,
        n_dnn_slices,
        1,
        nx,
        min_dnn_sum,
        min_bkg_neff,
        sig_mode,
        min_dnn_bkg_each,
        min_dnn_bkg_neff,
        min_bkg_frac,
    )
    dnn_slices = extend_outer_edges(dnn_slices, 0, nx + 1)

    result = []
    for xlo, xhi in dnn_slices:
        # ProjectionY attaches its result to gDirectory (here, the open output
        # file); detach so these die with the loop iteration instead of piling up
        # in the output file's in-memory object list.
        bkg_y_by_name = {
            name: [
                _detach(h.ProjectionY(f"_disc_{name}_{xlo}_{xhi}_{i}_y", xlo, xhi, "e"))
                for i, h in enumerate(hists)
            ]
            for name, hists in bkg2d_by_name.items()
        }
        sig_y = _detach(sig2d.ProjectionY(f"_disc_sig_{xlo}_{xhi}_y", xlo, xhi, "e"))
        n_bins = mass_bin_budget(bkg_y_by_name, 1, ny, max_mass_bins, bkg_per_mass_bin)
        mass_ranges = find_mass_bins(
            sig_y,
            bkg_y_by_name,
            n_bins,
            1,
            ny,
            min_mass_each,
            min_bkg_frac,
            min_mass_bkg_neff,
        )
        mass_ranges = extend_outer_edges(mass_ranges, 0, ny + 1)
        result.append({"x_range": (xlo, xhi), "y_ranges": mass_ranges})
    return result


def mass_bin_edges(y_axis, y_ranges):
    """Physical HME edges of the discovered y ranges, for booking the output TH1.

    The rebinned shapes used to be booked as n_bins over [0, n_bins], which threw
    the HME scale away and left every plot labelled by bin index. The ranges are
    contiguous and ordered, so the edges are each range's low edge plus the last
    range's upper edge. extend_outer_edges() pushes the outer ranges into the
    underflow/overflow bins, which have no finite edge of their own -- those are
    clamped to the axis limits.
    """
    n = y_axis.GetNbins()
    edges = [y_axis.GetBinLowEdge(max(lo, 1)) for lo, _ in y_ranges]
    edges.append(y_axis.GetBinUpEdge(min(y_ranges[-1][1], n)))
    return edges


def mkdir_titled(directory, path, title):
    """mkdir a nested path, putting `title` on the leaf directory only.

    TDirectory.mkdir() given a slashed path applies the title to the *first* level and
    hands the rest of the path down as the sub-levels' titles, so every parent ends up
    labelled with a stale fragment of whichever slice was created first -- the output
    currently has "muMu" titled "muMu/SR/res2b_dnn0". Walking the components keeps the
    parents clean and puts the label where it belongs.
    """
    parts = path.split("/")
    for part in parts[:-1]:
        directory = directory.GetDirectory(part) or directory.mkdir(part)
    return directory.mkdir(parts[-1], title)


def format_var_range(lo, hi, var):
    """One slice's edges on the sliced axis as a selection label, e.g. "1.20 < DNN < 4.50".

    `var` names that axis and comes from the configuration -- this script does not assume
    the analysis slices on a DNN score.

    Formatted here, by the code that discovered the edges, so the label travels with the
    histograms and nothing downstream has to re-derive it. Both edges open means the slice
    covers the whole axis, i.e. there is no selection to state -- that returns an empty
    string rather than a vacuous label.
    """
    if lo is None and hi is None:
        return ""
    if lo is None:
        return f"{var} < {hi:.2f}"
    if hi is None:
        return f"{var} > {lo:.2f}"
    return f"{lo:.2f} < {var} < {hi:.2f}"


def slice_ranges(x_axis, slices):
    """Physical edges of the discovered slices on the sliced axis, as [[lo, hi], ...].

    The slice x_ranges are bin indices, and extend_outer_edges() has already pushed the
    outermost ones into underflow/overflow -- those have no finite edge, so they are
    recorded as null and read back as an open-ended selection. Written out by run() so
    the plots can say which DNN selection each slice actually is; nothing downstream of
    the datacards needs it.
    """
    n = x_axis.GetNbins()
    ranges = []
    for sl in slices:
        lo, hi = sl["x_range"]
        ranges.append(
            [
                None if lo < 1 else x_axis.GetBinLowEdge(lo),
                None if hi > n else x_axis.GetBinUpEdge(hi),
            ]
        )
    return ranges


def rebin_hist_2d(hist2d, slices, name, naming):
    """Given the discovered slice structure, produce one final TH1 per slice
    for this specific histogram (nominal or a systematic variation)."""
    outputs = []
    for slice_idx, sl in enumerate(slices):
        xlo, xhi = sl["x_range"]
        edges = array.array("d", mass_bin_edges(hist2d.GetYaxis(), sl["y_ranges"]))
        # Detached for the same reason as the projections above: Write() targets
        # gDirectory regardless, so nothing needs these to stay attached.
        h = _detach(
            ROOT.TH1D(
                naming.name(name, slice_idx),
                name,
                len(edges) - 1,
                edges,
            )
        )
        for bin_idx, (ylo, yhi) in enumerate(sl["y_ranges"], start=1):
            err = array.array("d", [0.0])
            content = hist2d.IntegralAndError(xlo, xhi, ylo, yhi, err)
            h.SetBinContent(bin_idx, content)
            h.SetBinError(bin_idx, err[0])
        outputs.append(h)
    return outputs


def process_category(
    in_file,
    out_file,
    channel,
    category,
    cfg,
    mass,
    era,
    discovery_files,
    n_dnn_slices,
    max_mass_bins,
    min_dnn_sum,
    min_mass_each,
    min_signal,
    min_bkg_neff=0.0,
    min_bkg_frac=0.0,
    min_mass_bkg_neff=0.0,
    bkg_per_mass_bin=0.0,
    sig_mode="sb",
    min_dnn_bkg_each=0.0,
    min_dnn_bkg_neff=0.0,
):
    prefix = f"{channel}/{category}/"
    cat_dir = in_file.Get(f"{channel}/{category}")
    if not cat_dir:
        print(f"  [skip] {channel}/{category}: not found in {in_file.GetName()}")
        return
    keys = [k.GetName() for k in cat_dir.GetListOfKeys()]

    signal_keys = [
        applyParameters(pattern, {cfg["signal_param_name"]: mass})
        for pattern in cfg["signal_hist_name_patterns"]
    ]
    background_names = [
        hist_name
        for (base_name, hist_name, allowed_channels) in cfg["background_entries"]
        if not allowed_channels or channel in allowed_channels
    ]

    def load2d(f, key):
        return get_hist(f, prefix + key)

    disc_sig = sum_hists(
        [load2d(f, key) for f in discovery_files for key in signal_keys]
    )
    disc_bkg_by_name = {}
    for bkg_key in background_names:
        per_era = [load2d(f, bkg_key) for f in discovery_files]
        per_era = [h for h in per_era if h is not None]
        if len(per_era) == len(discovery_files):
            disc_bkg_by_name[bkg_key] = per_era

    sig_integral = disc_sig.Integral() if disc_sig is not None else 0
    if disc_sig is None or sig_integral < min_signal or not disc_bkg_by_name:
        if not disc_bkg_by_name:
            reason = "no background histograms found in all discovery eras"
        else:
            reason = f"signal too small for discovery ({sig_integral} < {min_signal})"
        print(f"    [skip] {channel}/{category} MX={mass}: {reason}, skipping")
        return

    slices = discover_binning(
        disc_sig,
        disc_bkg_by_name,
        n_dnn_slices,
        max_mass_bins,
        min_dnn_sum,
        min_mass_each,
        min_bkg_neff,
        min_bkg_frac,
        min_mass_bkg_neff,
        bkg_per_mass_bin,
        sig_mode,
        min_dnn_bkg_each,
        min_dnn_bkg_neff,
    )
    # The DNN selection each slice stands for is otherwise nowhere in the output: the
    # sliced categories are named by index and the surviving axis is the HME one. It is
    # the slice directory's own title, so it travels with the histograms and there is no
    # side-car path to hand a reader correctly and no key name for the two ends to agree
    # on -- anything that can open the shapes can already read it.
    naming = cfg["naming"]
    ranges = slice_ranges(disc_sig.GetXaxis(), slices)
    for slice_idx in range(len(slices)):
        mkdir_titled(
            out_file,
            f"{channel}/{naming.name(category, slice_idx)}",
            format_var_range(*ranges[slice_idx], var=cfg["slice_var"]),
        )

    for key in keys:
        hist2d = load2d(in_file, key)
        if hist2d is None or hist2d.GetDimension() != 2:
            continue
        rebinned = rebin_hist_2d(hist2d, slices, key, naming)
        for slice_idx, h in enumerate(rebinned):
            out_cat = naming.name(category, slice_idx)
            out_file.cd(f"{channel}/{out_cat}")
            h.Write(key)


def run(
    input_dir,
    output_dir,
    config_path,
    era,
    discovery_eras,
    n_dnn_slices,
    max_mass_bins,
    min_dnn_sum,
    min_mass_each,
    min_signal,
    min_bkg_neff=0.0,
    min_bkg_frac=0.0,
    min_mass_bkg_neff=0.0,
    bkg_per_mass_bin=0.0,
    sig_mode="sb",
    min_dnn_bkg_each=0.0,
    min_dnn_bkg_neff=0.0,
    category_pattern=None,
    slice_var="x",
):
    cfg = load_config(config_path)
    # Taken from the command line rather than from the configuration just read, so that an
    # override reaches the names actually written. HistRebinTask resolves the value once
    # and hands the same one to every reader, so the written names and the datacard bins
    # cannot come from different patterns.
    cfg["naming"] = CategoryNaming(category_pattern)
    cfg["slice_var"] = slice_var
    model = cfg["model"]

    for mass in cfg["mass_values"]:
        file_name = model.getInputFileName(era, {cfg["signal_param_name"]: mass})
        in_file = open_input_file(input_dir, model, era, mass, cfg["signal_param_name"])
        discovery_files = [
            open_input_file(input_dir, model, disc_era, mass, cfg["signal_param_name"])
            for disc_era in discovery_eras
        ]

        # output_dir is already period-scoped (HistRebinTask.output() = HistRebin/<era>),
        # so strip the leading "<era>/" that getInputFileName() adds -- otherwise the era
        # ends up doubled in the output path.
        era_prefix = era + "/"
        out_rel = (
            file_name[len(era_prefix) :]
            if file_name.startswith(era_prefix)
            else file_name
        )
        out_path = os.path.join(output_dir, out_rel)
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        out_file = ROOT.TFile.Open(out_path, "RECREATE")

        print(f"Rebinning era={era} MX={mass} -> {out_path}")
        for channel in cfg["channels"]:
            for category in cfg["categories"]:
                process_category(
                    in_file,
                    out_file,
                    channel,
                    category,
                    cfg,
                    mass,
                    era,
                    discovery_files,
                    n_dnn_slices,
                    max_mass_bins,
                    min_dnn_sum,
                    min_mass_each,
                    min_signal,
                    min_bkg_neff,
                    min_bkg_frac,
                    min_mass_bkg_neff,
                    bkg_per_mass_bin,
                    sig_mode,
                    min_dnn_bkg_each,
                    min_dnn_bkg_neff,
                )

        out_file.Close()
        in_file.Close()
        for f in discovery_files:
            f.Close()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Rebin 2D (DNN x HME) histograms into significance-sliced 1D shapes."
    )
    parser.add_argument(
        "--input",
        required=True,
        type=str,
        help="base directory containing <era>/<var>/<var>.root",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=str,
        help="output base directory, mirrors --input layout",
    )
    parser.add_argument(
        "--config", required=True, type=str, help="datacard configuration yaml"
    )
    parser.add_argument(
        "--era", required=True, type=str, help="era to rebin and write out"
    )
    parser.add_argument(
        "--discovery-eras",
        required=False,
        type=str,
        default=None,
        help="comma-separated eras summed to discover bin edges (defaults to --era alone)",
    )
    parser.add_argument(
        "--category-pattern",
        required=False,
        type=str,
        default=None,
        help="pattern naming the categories a base category is sliced into, e.g. "
        f"'{CategoryNaming.default_pattern}' (the default). Must use both "
        "{base_category} and {slice_idx}; every reader of these shapes is handed the "
        "same pattern, and it is what they parse the names back with",
    )
    parser.add_argument(
        "--slice-var",
        required=False,
        type=str,
        default="x",
        help="name of the sliced axis, used only to label each slice directory with the "
        "selection it stands for (e.g. 'DNN' -> '0.80 < DNN < 1.00')",
    )
    parser.add_argument(
        "--n-dnn-slices",
        required=False,
        type=int,
        default=4,
        help="fixed number of DNN slices per category (must be the same for every mass point)",
    )
    parser.add_argument(
        "--max-mass-bins",
        required=False,
        type=int,
        default=10,
        help="mass bins to aim for inside each DNN slice; backed off until achievable",
    )
    parser.add_argument(
        "--min-dnn-bkg-sum",
        required=False,
        type=float,
        default=1.0,
        help="minimum summed-background yield required in a DNN slice",
    )
    parser.add_argument(
        "--min-mass-bkg-each",
        required=False,
        type=float,
        default=0.01,
        help="minimum yield required of every individual background in a mass bin",
    )
    parser.add_argument(
        "--min-bkg-neff",
        required=False,
        type=float,
        default=0.0,
        help="minimum effective MC entries ((sum/err)^2) required of the summed "
        "background for a DNN slice boundary to be selectable",
    )
    parser.add_argument(
        "--min-bkg-frac",
        required=False,
        type=float,
        default=0.0,
        help="backgrounds contributing less than this fraction of the total are "
        "exempt from --min-mass-bkg-each, so a statistically starved minor "
        "process cannot collapse a slice to a single mass bin",
    )
    parser.add_argument(
        "--min-signal",
        required=False,
        type=float,
        default=0.5,
        help="minimum discovery signal yield required to slice a category at all "
        "(below this, e.g. boosted at low MX, the category is skipped entirely)",
    )
    parser.add_argument(
        "--min-mass-bkg-neff",
        required=False,
        type=float,
        default=0.0,
        help="minimum effective MC entries required of the summed background in "
        "every mass bin. The --min-bkg-neff analogue for the bins inside a "
        "slice, which is where essentially every fit bin lives",
    )
    parser.add_argument(
        "--bkg-per-mass-bin",
        required=False,
        type=float,
        default=0.0,
        help="target summed-background yield per mass bin; caps the bin count at "
        "B_slice/this instead of always using --max-mass-bins. 0 disables, "
        "restoring the fixed --max-mass-bins for every slice",
    )
    parser.add_argument(
        "--min-dnn-bkg-each",
        required=False,
        type=float,
        default=0.0,
        help="minimum yield required of every non-negligible background in a DNN "
        "slice. The --min-mass-bkg-each analogue for the slice boundaries "
        "themselves, which were previously gated on the summed background only: "
        "a background that had fluctuated negative was hidden inside a healthy "
        "total, and the resulting slice could not be made into a datacard. "
        "Backgrounds below --min-bkg-frac of the category total are exempt. "
        "0 disables",
    )
    parser.add_argument(
        "--min-dnn-bkg-neff",
        required=False,
        type=float,
        default=0.0,
        help="minimum effective MC entries required of every non-negligible "
        "background in a DNN slice (same exemption as --min-dnn-bkg-each). Much "
        "stronger than --min-dnn-bkg-each and correspondingly expensive: on "
        "Run3_Early a threshold of 4 cost 35% of the combined Asimov Z in "
        "muMu/SR/res2b at m500, against 2.3% for the yield floor alone. 0 disables",
    )
    parser.add_argument(
        "--significance-mode",
        required=False,
        type=str,
        default="sb",
        choices=["sb", "asimov"],
        help="figure of merit for DNN slice boundaries: 'sb' = S/sqrt(B+sigmaB^2), "
        "'asimov' = Poisson-correct Asimov significance (valid at low B)",
    )
    args = parser.parse_args()

    discovery_eras = (
        args.discovery_eras.split(",") if args.discovery_eras else [args.era]
    )

    run(
        args.input,
        args.output,
        args.config,
        args.era,
        discovery_eras,
        args.n_dnn_slices,
        args.max_mass_bins,
        args.min_dnn_bkg_sum,
        args.min_mass_bkg_each,
        args.min_signal,
        args.min_bkg_neff,
        args.min_bkg_frac,
        args.min_mass_bkg_neff,
        args.bkg_per_mass_bin,
        args.significance_mode,
        args.min_dnn_bkg_each,
        args.min_dnn_bkg_neff,
        category_pattern=args.category_pattern,
        slice_var=args.slice_var,
    )
