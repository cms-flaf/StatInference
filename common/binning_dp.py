"""Prefix-sum tables over a 2D shape, so the window search can score any bin in O(1).

The window search tries ~ny^2 windows, each with an exact partition search over ~nx^2
candidate bins, which is only affordable with table lookups instead of TH2 integrals.
The search is written with the window on y; build_cells transposes the input for a
window on x.

The gate (_bin_passes) and the figure of merit (significance) stay in rebin_2d.py and are
passed in, so each has one definition.

Notes:
- Under/overflow are included (arrays are indexed 0..n+1), to match what the TH2
  integrals in rebin_2d.py return.
- Differencing prefix sums can turn an exact 0 into -2e-16, which the `value < 0`
  positivity gate would reject. Yields are snapped to 0 below an epsilon scaled by each
  process's L1 norm (~1e-8 events, far below the 0.01 floor).
- check_window_mask() re-checks a sample of every chosen window against _bin_passes and
  significance on each run.
"""

import math

import numpy as np


def _hist_arrays(hist):
    """(values, variances) as (nx+2, ny+2) arrays indexed [binx, biny], bins 0..n+1."""
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

    Scaled by the L1 norm, not the sum: with negative weights a process can nearly cancel
    to a small sum while its prefix sums still lose precision on the L1 scale.
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

    Backgrounds are summed across discovery eras (errors in quadrature), as
    _bkg_yields()/_bkg_errors() do.
    """

    def __init__(self, nx, ny, sig, bkg, var, significance):
        self.nx = nx
        self._significance = significance
        self.ny = ny
        self.names = list(bkg)
        self._sig = _prefix(sig)
        self._bkg = {name: _prefix(bkg[name]) for name in self.names}
        self._var = {name: _prefix(var[name]) for name in self.names}
        # no backgrounds at all is possible (a replay whose input lacks them); keep arrays
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
        """significance()^2 for one rectangle; Z^2 adds across bins, so a partition's
        value is the sum of its bins'."""
        s = self.signal(x0, x1, y0, y1)
        b = self.total_bkg(x0, x1, y0, y1)
        b_err = self.total_bkg_error(x0, x1, y0, y1)
        return self._significance(s, b, b_err, mode) ** 2

    def exempt(self, x0, x1, y0, y1, min_frac):
        """minor_backgrounds() over the rectangle about to be subdivided."""
        if min_frac <= 0:
            return set()
        y = self.yields(x0, x1, y0, y1)
        total = sum(y.values())
        if total <= 0:
            return set()
        return {n for n, v in y.items() if v < min_frac * total}


def build_cells(sig2d, bkg2d_by_name, significance, transpose=False):
    """Cells for one channel/category.

    sig2d is the summed discovery signal; bkg2d_by_name is {background: [hist per era]}.
    transpose swaps x and y, which is how the window search runs with the window on x.
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

    Returns (values, partitions): values[k] is the best k-part total score (-inf if none)
    and partitions[k] its (lo, hi) index pairs; index 0 is unused. Every k is returned
    because the bin count is chosen from the whole curve.

    Exact interval DP, O(max_parts * n^2), valid because the score is additive over parts.
    Invalid ranges are excluded outright, so the background gates are hard constraints.
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
    """(total Z^2, per-slice Z^2) of a finished binning, from the recorded ranges.

    Computed for every strategy, so binnings can be compared from binning.json alone.
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
    """M[i, j] = the sum over bins (i+1)..(j+1), from one column of a prefix sum.

    The whole (n, n) table of candidate ranges is one outer subtraction. With
    include_outer the first and last bins also take the under/overflow.
    """
    lo = column[1 : n + 1].copy()
    hi = column[2 : n + 2].copy()
    if include_outer:
        lo[0] = column[0]
        hi[-1] = column[n + 2]
    return hi[None, :] - lo[:, None]


def _neff_matrix(value, error):
    """effective_entries() over arrays. error <= 0 gives inf if value > 0, else 0, as
    the scalar does (a nan would silently fail every threshold)."""
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
    """Give back every bin whose split raised Z^2 by less than `threshold`.

    Z^2 almost always grows with the bin count, so without this the search would use
    every bin it is allowed. The caller sets `threshold` as a fraction of the category's
    total, and the count is stepped down until the last bin pays for itself.
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
    """(score, valid) for every candidate x bin inside the window y0..y1.

    _bin_passes() and significance() restated as array operations: calling them per
    candidate is far too slow for the search. check_window_mask() guards against the two
    versions drifting apart.

    `exempt` (processes excused from the per-bin floors) is judged once over the whole
    window. Positivity is never excused.
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
    """Re-derive a random sample of window_tables() with _bin_passes and significance(),
    and stop the run on any disagreement."""
    nx = cells.nx
    rng = np.random.default_rng(seed)
    for _ in range(n_samples):
        a = int(rng.integers(1, nx + 1))
        b = int(rng.integers(a, nx + 1))
        # the outermost bins include the x under/overflow, as in window_tables
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
