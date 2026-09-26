"""Prefix-sum cells over a 2D shape, so any candidate binning can be scored in O(1).

The greedy binner asks its questions one candidate at a time and re-integrates the
histograms for each, which is affordable because it only ever looks at ~nx windows. The
window search looks at ~ny^2 windows and runs an exact partition search over ~nx^2
candidate bins of x inside every one (in the cells' own frame, where the window is on y;
build_cells transposes the input when the window is on x), so the integrals have to become table lookups or
the search is not worth running.
That is all this module does: it turns the TH2s into cumulative arrays and hands back the
same yields, errors and figure of merit the greedy path would have computed.

It deliberately knows nothing about what a valid bin is or how a bin is scored. The
gate (_bin_passes) and the figure of merit (significance) live in the 2D binner and are
passed in as callables, so there is exactly one definition of each and this module
cannot drift from it. What it provides is the *arguments* those gates take -- the
yields and errors dicts -- built from the tables instead of from ROOT.

Under- and overflow are carried, and that is load-bearing rather than tidy. _integral()
in the 2D binner reads a TH2 with Integral(lo, hi, 0, -1), and ROOT reads biny2 < biny1 as
"the whole y axis including under- and overflow", so every slice-level quantity in the
greedy path already includes the y under/overflow rows. A table that stopped at bin ny
would disagree with the greedy path by a fraction of an event, which is easily enough to
move a boundary and produce a difference nobody can explain. So the arrays are indexed
0..n+1 throughout, and a slice-level query spans y = 0..ny+1.

Differencing prefix sums is not exact, and one of the ways it is inexact has teeth.
A window whose content is ~1e8 times smaller than the whole plane's -- routine for DY in
a narrow window at high y -- loses about eight digits to cancellation, so a quantity that
is exactly 0.0 in ROOT comes back as -2e-16 here. That is numerically nothing and
physically nothing, but _bin_passes() tests `value < 0` to reject a background that has
gone negative, and on the real inputs it rejected 137 sampled windows on the strength
of a rounding sign. So every yield is snapped to zero below a per-process epsilon scaled
by that process's own L1 norm, which is the size of the cancellation that can occur. The
epsilon lands around 1e-8 events for the largest process here, against a magnitude floor
(min_bin_bkg_each) of 0.01 -- six orders of margin, so the snap cannot mask a real
negative.

The check that all of this is true is check_window_mask(), not inspection: every
window the search settles on has a sample of its table entries re-derived with _bin_passes and
significance on the histograms themselves, and the run stops on any disagreement.
"""

import math

import numpy as np


def _hist_arrays(hist):
    """(values, variances) as (nx+2, ny+2) arrays indexed [binx, biny], bins 0..n+1.

    Read through GetBinContent/GetBinError rather than the internal buffer: the buffer
    is faster but its layout is a ROOT implementation detail, and getting it wrong is
    silent. This runs once per histogram per category, which is nothing next to the
    search it feeds.
    """
    nx = hist.GetNbinsX()
    ny = hist.GetNbinsY()
    values = np.empty((nx + 2, ny + 2), dtype=np.float64)
    variances = np.empty((nx + 2, ny + 2), dtype=np.float64)
    for bx in range(nx + 2):
        for by in range(ny + 2):
            values[bx, by] = hist.GetBinContent(bx, by)
            err = hist.GetBinError(bx, by)
            variances[bx, by] = err * err
    return values, variances


def _prefix(a):
    """Inclusive 2D prefix sums, shape (nx+3, ny+3), with P[i, j] = sum(a[:i, :j])."""
    p = np.zeros((a.shape[0] + 1, a.shape[1] + 1), dtype=np.float64)
    p[1:, 1:] = a.cumsum(axis=0).cumsum(axis=1)
    return p


def _cancellation_eps(a, rel=1e-12):
    """How large a value this array's prefix sums can invent out of rounding.

    Scaled by the L1 norm rather than the sum, because the cancellation that matters
    comes from negative-weight events: a process whose bins nearly cancel has a small
    sum and a large L1, and it is the L1 that sets how much precision the difference of
    two prefix sums can lose. rel is a few hundred times the double-precision epsilon,
    which covers accumulating over the ~1e4 cells of one of these planes.
    """
    return rel * float(np.abs(a).sum())


def _snap(value, eps):
    """Zero out a value that is indistinguishable from zero at this array's precision."""
    return 0.0 if abs(value) < eps else value


def _rect(p, x0, x1, y0, y1):
    """Sum over the inclusive bin-index rectangle [x0..x1] x [y0..y1]."""
    return float(p[x1 + 1, y1 + 1] - p[x0, y1 + 1] - p[x1 + 1, y0] + p[x0, y0])


class Cells:
    """Every quantity the binning gates and the figure of merit need, as O(1) lookups.

    Background yields are summed across discovery eras and their MC errors added in
    quadrature, which is what _bkg_yields()/_bkg_errors() do -- the eras are summed
    because the binning is derived from the combination it will be applied to.
    """

    def __init__(self, nx, ny, sig, bkg, var, significance):
        self.nx = nx
        self._significance = significance
        self.ny = ny
        self.names = list(bkg)
        self._sig = _prefix(sig)
        self._bkg = {name: _prefix(bkg[name]) for name in self.names}
        self._var = {name: _prefix(var[name]) for name in self.names}
        # An empty background dict is reachable -- a replay whose input is missing a
        # process it was derived with -- and sum() of nothing is the integer 0, which is
        # not an array. Zeros keep every query defined and answering zero, which is the
        # truthful answer when there is no background to integrate.
        zero = np.zeros_like(sig)
        bkg_tot = sum(bkg.values()) if bkg else zero
        var_tot = sum(var.values()) if var else zero
        self._bkg_tot = _prefix(bkg_tot)
        self._var_tot = _prefix(var_tot)
        self._eps_sig = _cancellation_eps(sig)
        self._eps_bkg = {name: _cancellation_eps(bkg[name]) for name in self.names}
        self._eps_tot = _cancellation_eps(bkg_tot)
        self._eps_var = {name: _cancellation_eps(var[name]) for name in self.names}
        self._eps_var_tot = _cancellation_eps(var_tot)

    def signal(self, x0, x1, y0, y1):
        return _snap(_rect(self._sig, x0, x1, y0, y1), self._eps_sig)

    def total_bkg(self, x0, x1, y0, y1):
        return _snap(_rect(self._bkg_tot, x0, x1, y0, y1), self._eps_tot)

    def total_bkg_error(self, x0, x1, y0, y1):
        return math.sqrt(max(_rect(self._var_tot, x0, x1, y0, y1), 0.0))

    def yields(self, x0, x1, y0, y1):
        """{background: yield}, the dict the gates take."""
        return {
            n: _snap(_rect(self._bkg[n], x0, x1, y0, y1), self._eps_bkg[n])
            for n in self.names
        }

    def errors(self, x0, x1, y0, y1):
        """{background: MC statistical error}, the dict the gates take."""
        return {
            n: math.sqrt(max(_rect(self._var[n], x0, x1, y0, y1), 0.0))
            for n in self.names
        }

    def score(self, x0, x1, y0, y1, mode):
        """The figure of merit for one cell, squared.

        Squared because Asimov Z^2 is what adds across independent counting bins, so a
        partition's value is the sum of its cells' scores and the search can be a
        dynamic program. `significance` is the 2D binner's own, passed in, so the cells
        are ranked by exactly the quantity the greedy path maximises.
        """
        s = self.signal(x0, x1, y0, y1)
        b = self.total_bkg(x0, x1, y0, y1)
        b_err = self.total_bkg_error(x0, x1, y0, y1)
        return self._significance(s, b, b_err, mode) ** 2

    def exempt(self, x0, x1, y0, y1, min_frac):
        """minor_backgrounds() over a rectangle, for the range the caller is about to
        subdivide -- never for the candidate sub-range itself, which is circular.
        """
        if min_frac <= 0:
            return set()
        y = self.yields(x0, x1, y0, y1)
        total = sum(y.values())
        if total <= 0:
            return set()
        return {n for n, v in y.items() if v < min_frac * total}


def build_cells(sig2d, bkg2d_by_name, significance, transpose=False):
    """Cells for one channel/category, from the same histograms discover_binning() reads.

    sig2d is the already-summed discovery signal; bkg2d_by_name is
    {background: [hist per discovery era]}, summed here the way _bkg_yields() sums it.

    transpose swaps the two axes, so the cells' x is the input's y and vice versa. The
    window search is written with the window on y; this is how it runs on x instead.
    """
    nx = sig2d.GetNbinsX()
    ny = sig2d.GetNbinsY()
    sig, _ = _hist_arrays(sig2d)
    bkg = {}
    var = {}
    for name, hists in bkg2d_by_name.items():
        values = np.zeros((nx + 2, ny + 2), dtype=np.float64)
        variances = np.zeros((nx + 2, ny + 2), dtype=np.float64)
        for h in hists:
            v, e2 = _hist_arrays(h)
            values += v
            variances += e2
        bkg[name] = values
        var[name] = variances
    if transpose:
        nx, ny = ny, nx
        sig = sig.T.copy()
        bkg = {name: a.T.copy() for name, a in bkg.items()}
        var = {name: a.T.copy() for name, a in var.items()}
    return Cells(nx, ny, sig, bkg, var, significance)


NEG = -np.inf


def partition_dp(score, valid, max_parts):
    """Best partition of [0..n-1] into exactly k contiguous valid ranges, for every k.

    Returns (values, partitions): values[k] is the total score of the best k-part
    partition, or -inf if there is none, and partitions[k] is that partition as a list of
    (lo, hi) index pairs. Index 0 is unused in both so that k reads as the part count.

    Every k is returned rather than just the requested one, because the caller needs the
    whole curve and it costs nothing to keep: the bin count is chosen by comparing what
    the partition scores at each count, and trim_by_marginal_gain() gives a bin back by
    comparing values[k] against values[k-1]. Computing them one at a time would re-run the
    same DP.

    Exact, not greedy. The score is additive over the parts -- Asimov Z^2 is what adds
    across independent counting bins -- so the optimal k-part partition of a prefix is
    built from an optimal (k-1)-part partition of a shorter prefix, and the usual
    interval DP applies. O(max_parts * n^2).

    `valid` is a mask, not a penalty: an invalid range is unreachable rather than merely
    expensive, so a partition containing one cannot be returned at any score. That is what
    keeps the background gates hard constraints instead of preferences.
    """
    n = score.shape[0]
    max_parts = max(1, min(max_parts, n))
    masked = np.where(valid, score, NEG)

    dp = np.full((max_parts + 1, n + 1), NEG)
    arg = np.full((max_parts + 1, n + 1), -1, dtype=int)
    dp[0, 0] = 0.0
    for k in range(1, max_parts + 1):
        for e in range(k, n + 1):
            # candidate: the last part is [s .. e-1], the rest is a (k-1)-part prefix
            cand = dp[k - 1, :e] + masked[:e, e - 1]
            s = int(np.argmax(cand))
            if cand[s] > NEG:
                dp[k, e] = cand[s]
                arg[k, e] = s

    values = [NEG] * (max_parts + 1)
    partitions = [None] * (max_parts + 1)
    for k in range(1, max_parts + 1):
        if dp[k, n] == NEG:
            continue
        values[k] = float(dp[k, n])
        parts = []
        e = n
        for kk in range(k, 0, -1):
            s = int(arg[kk, e])
            parts.append((s, e - 1))
            e = s
        partitions[k] = list(reversed(parts))
    return values, partitions


def binning_objective(cells, slices, mode):
    """(total Z^2, per-slice Z^2) for a finished binning.

    The quantity the search maximises, recomputed from the binning that was actually
    written rather than carried out of the optimiser. That makes it meaningful for the
    greedy strategy too -- which never computes it -- so the two can be compared from
    their binning.json alone, and it means a discrepancy between what the optimiser
    thought it achieved and what the shapes contain shows up as a discrepancy rather
    than going unnoticed.

    Note this is evaluated on the ranges as recorded, i.e. after extend_outer_edges has
    pushed the outermost ones into under/overflow, so it describes the shapes on disk.
    """
    per_slice = []
    for sl in slices:
        if "y_range" in sl:  # selection on y, bins along x
            ylo, yhi = sl["y_range"]
            per_slice.append(
                sum(cells.score(a, b, ylo, yhi, mode) for a, b in sl["x_ranges"])
            )
        else:  # selection on x, bins along y
            xlo, xhi = sl["x_range"]
            per_slice.append(
                sum(cells.score(xlo, xhi, a, b, mode) for a, b in sl["y_ranges"])
            )
    return sum(per_slice), per_slice


def _ranges_from_column(column, n, include_outer=False):
    """M[i, j] = the rectangle sum for bins (i+1)..(j+1), from one column of a prefix sum.

    For a fixed x window, P[x1+1, :] - P[x0, :] is the running sum along y, and every
    y range is a difference of two of its entries -- so the whole (n, n) table of
    candidate ranges is one outer subtraction rather than n^2 lookups.

    With include_outer the first and last bins swallow the under- and overflow, which is
    what extend_outer_edges() does to the outermost ranges before they are written. A
    search that scores the ranges without them is optimising something slightly different
    from what ends up in the datacard -- 0.11% of the background on the DNN axis of the
    HH->bbWW shapes, small but not nothing, and concentrated entirely in the two
    outermost slices.
    """
    lo = column[1 : n + 1].copy()
    hi = column[2 : n + 2].copy()
    if include_outer:
        lo[0] = column[0]
        hi[-1] = column[n + 2]
    return hi[None, :] - lo[:, None]


def _neff_matrix(value, error):
    """effective_entries() over arrays, with its two degenerate cases kept.

    error <= 0 means the yield carries no MC uncertainty at all: infinite effective
    entries if there is something there, none if there is not. Taking (v/e)^2 blindly
    would make the first case a division by zero and the second a nan, and a nan
    compares false against every threshold -- which would silently reject exactly the
    empty bins the scalar gate accepts.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(error > 0, value / np.where(error > 0, error, 1.0), 0.0) ** 2
    degenerate = np.where(value > 0, np.inf, 0.0)
    return np.where(error > 0, ratio, degenerate)


def _asimov_z2_matrix(s, b, b_err):
    """asimov_significance(...)**2 over arrays, matching the scalar branch for branch."""
    out = np.zeros_like(s)
    live = (b > 0) & (s > 0)
    if not live.any():
        return out
    var = b_err**2
    no_err = live & (b_err <= 0)
    if no_err.any():
        ss, bb = s[no_err], b[no_err]
        out[no_err] = np.maximum(2.0 * ((ss + bb) * np.log1p(ss / bb) - ss), 0.0)
    with_err = live & (b_err > 0)
    if with_err.any():
        ss, bb, vv = s[with_err], b[with_err], var[with_err]
        term1 = (ss + bb) * np.log(((ss + bb) * (bb + vv)) / (bb * bb + (ss + bb) * vv))
        term2 = (bb * bb / vv) * np.log1p(vv * ss / (bb * (bb + vv)))
        out[with_err] = np.maximum(2.0 * (term1 - term2), 0.0)
    return out


def _sb_z2_matrix(s, b, b_err):
    """significance(..., mode="sb")**2 over arrays: S^2 / (B + sigma_B^2)."""
    denom = b + b_err**2
    ok = (b > 0) & (denom > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(ok, s**2 / np.where(ok, denom, 1.0), 0.0)


def trim_by_marginal_gain(values_per_slice, counts, threshold):
    """Give back every bin whose last split bought less than `threshold`.

    With Z^2 very nearly monotone in bin count, a search told only to maximise it will
    always spend whatever it is given. So the binding constraint has to be the value of a
    bin, not the size of a pot. A bin is kept only if the split that created it raised
    this category's Z^2 by at least `threshold` -- expressed as a fraction of the
    category's own achievable total, so it means "worth something on the scale that
    matters here" rather than a yield in events. A bin that buys nothing is
    MC-statistical exposure bought for nothing; this is what declines it.

    Stepping down one bin at a time and stopping at the first marginal gain that clears
    the threshold is exact for a concave curve, which Z^2 in bin count is: each extra
    split has less left to separate. It is not assumed -- a curve that is not concave
    simply stops at its first qualifying step, which is still a count whose last bin paid
    for itself.
    """
    if threshold <= 0:
        return list(counts)
    trimmed = []
    for values, n in zip(values_per_slice, counts):
        while n > 1 and values[n] - values[n - 1] < threshold:
            n -= 1
        trimmed.append(n)
    return trimmed


def _axis_ranges(prefix, a, b, n, axis, include_outer=False):
    """All ranges along `axis`, at a fixed window [a, b] on the other one."""
    if axis == "y":
        column = prefix[b + 1] - prefix[a]
    else:
        column = prefix[:, b + 1] - prefix[:, a]
    return _ranges_from_column(column, n, include_outer)


def window_tables(cells, y0, y1, exempt, knobs, mode):
    """(score, valid) over every candidate bin of the *sliced-on-y* layout.

    The window is a cut on the cells' y around the signal, and the bins run along their
    x (build_cells transposes the input when the window is on x). The gates are
    _bin_passes() and the figure of merit significance(), as
    whole-array operations: the scalar path -- building a yields dict and an errors dict
    per candidate and calling _bin_passes -- costs about 80 microseconds a candidate, and
    there are nx(nx+1)/2 of them per window, for every window the search tries.

    This is the one place that restates the gate rather than calling it, which is a real
    risk: two spellings of the same rule are two rules as soon as one of them is edited.
    It is contained by check_window_mask(), which re-evaluates the canonical _bin_passes on a
    random sample of candidates and raises on the first disagreement -- so the fast path
    is checked against the slow one on real inputs, on every run, rather than at review
    time.

    The bins tile the whole x axis including its under/overflow, because everything inside
    the window is kept and binned. The window itself is a selection: what falls outside it is
    discarded, not swept into an outer bin, which is the difference between a cut and a
    slice.

    `exempt` is the set of processes excused from the per-bin floors, judged once over
    the whole window -- see minor_backgrounds(). Positivity is never excused -- see
    _bin_passes.
    """
    nx = cells.nx
    ranges = lambda prefix: _axis_ranges(prefix, y0, y1, nx, "x", True)

    sig = ranges(cells._sig)
    b_tot = ranges(cells._bkg_tot)
    b_tot = np.where(np.abs(b_tot) < cells._eps_tot, 0.0, b_tot)
    v_tot = ranges(cells._var_tot)
    err = np.sqrt(np.maximum(v_tot, 0.0))

    valid = np.triu(np.ones((nx, nx), dtype=bool))
    if knobs["min_bin_bkg_neff"] > 0:
        valid &= _neff_matrix(b_tot, err) >= knobs["min_bin_bkg_neff"]
    for name in cells.names:
        b_p = ranges(cells._bkg[name])
        b_p = np.where(np.abs(b_p) < cells._eps_bkg[name], 0.0, b_p)
        valid &= b_p >= 0
        if name in exempt:
            continue
        magnitude = b_p > knobs["min_bin_bkg_each"]
        required = knobs["min_bin_bkg_each_neff"]
        if required > 0:
            v_p = ranges(cells._var[name])
            magnitude &= _neff_matrix(b_p, np.sqrt(np.maximum(v_p, 0.0))) >= required
        valid &= magnitude
    score = (
        _asimov_z2_matrix(sig, b_tot, err)
        if mode == "asimov"
        else _sb_z2_matrix(sig, b_tot, err)
    )
    return np.where(valid, score, 0.0), valid


def check_window_mask(
    cells,
    y0,
    y1,
    valid,
    score,
    exempt,
    bin_passes,
    mode,
    n_samples=100,
    seed=0,
):
    """The cross-check of window_tables() against the canonical scalar gate.

    Cheap -- a hundred calls against the thousands the tables replace -- and it is what
    makes the fast path safe to trust: the tables are re-derived on a random sample with
    _bin_passes and significance(), so if the two spellings of the rule ever disagree
    about a candidate bin the run stops here instead of quietly producing a different
    binning.
    """
    nx = cells.nx
    rng = np.random.default_rng(seed)
    for _ in range(n_samples):
        a = int(rng.integers(1, nx + 1))
        b = int(rng.integers(a, nx + 1))
        # the outermost bins carry the x under/overflow, because window_tables builds them
        # with include_outer -- the reference has to ask the same question
        lo = 0 if a == 1 else a
        hi = nx + 1 if b == nx else b
        want = bool(bin_passes(lo, hi, y0, y1, exempt))
        got = bool(valid[a - 1, b - 1])
        if got != want:
            raise AssertionError(
                f"binning_dp window gate disagrees with _bin_passes on bins {a}..{b} "
                f"of the window {y0}..{y1}: fast path says {'valid' if got else 'invalid'}, "
                f"_bin_passes says {'valid' if want else 'invalid'}."
            )
        if want:
            expected = cells.score(lo, hi, y0, y1, mode)
            if abs(expected - score[a - 1, b - 1]) > 1e-9 * max(abs(expected), 1.0):
                raise AssertionError(
                    f"binning_dp window score disagrees with significance() on bins "
                    f"{a}..{b} of the window {y0}..{y1}: {score[a - 1, b - 1]!r} vs {expected!r}."
                )
