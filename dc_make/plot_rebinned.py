"""Stacked plots of the rebinned 1D shapes that actually enter the datacards.

HistPlotTask plots the pre-rebinning histograms straight out of HistMergerTask, one
era at a time. That is the wrong object to look at once HistRebinTask has sliced the
2D plane into per-slice 1D shapes: what the fit sees is the
*rebinned* shape, summed over the meta-era's sub-eras exactly as maker.py's
getCombinedShape does at datacard-build time.

This script reproduces HistPlotTask's look (same PlotKit Plotter, same
config/plot/*.yaml style files, same signal-overlay convention) but takes the
rebinned files as input and sums the sub-eras first, so every panel corresponds
one-to-one with a datacard bin block.

One file per base category and mass, laid out as a grid: slices across,
lepton channels down. The fit sees these bins together, so they are easier to judge
together -- a slice that looks reasonable in eMu and pathological in eE is obvious
side by side and invisible in separate files.
"""

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

from StatInference.common.tools import (
    importROOT,
    CategoryNaming,
)
from StatInference.common.param_parse import extractParameters, applyParameters
from StatInference.dc_make.model import Model

ROOT = importROOT()

# Fallback colours, used only when a process is missing from processes.yaml so that
# a config gap degrades to an ugly plot rather than a crash.
FALLBACK_COLORS = {
    "TT": "kAzure+1",
    "DY": "kOrange-3",
    "ST": "kGreen+1",
    "VV": "kViolet-4",
    "W": "kRed-7",
    "SingleHiggs": "kGray+1",
}
FALLBACK_SIGNAL_COLOR = "kRed"


def load_config(config_path, n_slices=None, category_pattern=None):
    """Like hist_rebin_2d.load_config, but expands the base categories into the
    *sliced* names -- those are the directories that exist in the rebinned files and
    the ones that map onto datacard bins."""
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    # HistRebinTask's resolved value is authoritative when the caller passes it; a
    # command-line override would otherwise leave this deriving a different slice count
    # than the rebinned files contain. No slice count anywhere (or an explicit 0) means
    # the input is already binned and its categories are the ones listed -- same rule as
    # DatacardMaker, so the panels and the datacard bins cannot disagree.
    n_slices = n_slices
    if n_slices is None:
        n_slices = (cfg.get("binning") or {}).get("n_slices")
    # Same rule for the pattern naming those slices: the caller's resolved value wins,
    # otherwise the configuration's own.
    naming = (
        CategoryNaming(category_pattern)
        if category_pattern
        else CategoryNaming.fromConfig(cfg)
    )

    # multiple signal processes are allowed; each is drawn as its own stack entry
    signal_hist_names = []
    mass_values = None
    background_entries = []
    for entry in cfg["processes"]:
        if type(entry) == str:
            background_entries.append((entry, []))
            continue
        if entry.get("is_data", False):
            continue
        hist_name = entry.get("hist_name", entry["process"])
        if entry.get("is_signal", False):
            signal_hist_names.append(hist_name)
            if mass_values is None:
                mass_values = entry["param_values"]
        else:
            background_entries.append((hist_name, entry.get("channels", [])))

    if not signal_hist_names:
        raise RuntimeError("No signal process found in config")

    return {
        "model": Model.fromConfig(cfg["model"]),
        "channels": cfg["channels"],
        "naming": naming,
        "categories": (
            naming.expand(cfg["categories"], int(n_slices))
            if n_slices
            else list(cfg["categories"])
        ),
        "signal_hist_name_patterns": signal_hist_names,
        "signal_param_name": extractParameters(signal_hist_names[0])[0],
        "mass_values": mass_values,
        "background_entries": background_entries,
    }


def get_process_styles(ana_path, period, version, process_names):
    """Plot label + colour per process, taken from the analysis processes.yaml so
    the rebinned plots use the same colours as HistPlotTask."""
    styles = {}
    try:
        import FLAF.Common.Setup as Setup

        setup = Setup.Setup(ana_path, period, version)
        parent = setup.parent_processes
        for name in process_names:
            if name in parent:
                proc = parent[name]
                styles[name] = (proc.get("name", name), proc["color"])
    except Exception as e:
        print(
            f"Warning: could not load process styles from Setup ({e}); using fallbacks"
        )

    for name in process_names:
        if name not in styles or styles[name][1] == "kBlack":
            # kBlack is what Setup assigns to meta-process members it was not asked
            # to plot (every signal mass but the three in `to_plot`), which would be
            # invisible against the axis -- give the overlaid signal its own colour.
            is_signal = name not in FALLBACK_COLORS
            color = FALLBACK_SIGNAL_COLOR if is_signal else FALLBACK_COLORS[name]
            label = styles[name][0] if name in styles else name
            styles[name] = (label, color)
    return styles


def get_rebinned_axis_title(ana_path, period, version, variable):
    """Axis title for the surviving 1D axis, taken from the analysis's own plot config.

    The rebinned shapes are the y projection of the 2D `variable`, and histograms.yaml
    already describes both of its axes: the 2D entry's ``var_list`` is [x, y], and the y
    variable's own entry carries the ``x_title`` HistPlotTask draws it with. Reading it
    from there keeps this script free of any one analysis's axis label and keeps the
    rebinned plots labelled exactly like the pre-rebinning ones.

    Falls back to the variable's own name rather than a guess -- a plot with a plain
    label is a cosmetic problem, and inventing one would be a wrong one.
    """
    try:
        import FLAF.Common.Setup as Setup

        hists = Setup.Setup(ana_path, period, version).hists
        y_var = hists[variable]["var_list"][1]
        return hists[y_var]["x_title"]
    except Exception as e:
        print(f"Warning: no axis title for {variable} ({e}); using the variable name")
        return variable


def build_hist_cfg(x_title, y_title, log_y):
    """PlotKit reads binning/axis metadata out of histograms.yaml keyed by variable
    name. The rebinned axis is a bin index with no config entry (and needs none --
    the binning is already final), so hand it a synthetic one."""
    return {
        "rebinned": {
            "x_title": x_title,
            "y_title": y_title,
            "use_log_y": log_y,
            # Enough headroom for the legend without leaving the stack squashed into
            # the bottom decade (HistPlotter's 2000x is tuned for a full-page canvas).
            "max_y_sf": 100.0 if log_y else 1.5,
            "divide_by_bin_width": False,
        }
    }


class PanelBackend:
    """Draws a PlotKit StackSpec into an axes we already own.

    PlotKit's MplhepBackend makes its own figure in _new_figure() and saves it in
    render_stacked(). Swapping just those two things out lets every slice share
    one canvas while the drawing code -- and therefore the styling -- stays exactly
    the one HistPlotTask uses.
    """

    def __init__(self, ax, draw_cms_label, draw_legend):
        from FLAF.PlotKit.backends import MplhepBackend

        self._ax = ax
        self._draw_cms_label = draw_cms_label
        self._draw_legend = draw_legend

        outer = self

        class _Panel(MplhepBackend):
            def _new_figure(self, spec):
                import mplhep as hep

                class _NoSave:
                    def savefig(self, *args, **kwargs):
                        pass

                class _NoClose:
                    @staticmethod
                    def close(*args, **kwargs):
                        pass

                return _NoClose(), hep, _NoSave(), outer._ax, None

            def _cms_label(self, hep, ax, spec):
                # Only the leftmost panel carries "CMS Simulation" + the lumi.
                if outer._draw_cms_label:
                    super()._cms_label(hep, ax, spec)

        self._impl = _Panel()

    def render(self, spec):
        self._impl.render_stacked(spec, None)
        if not self._draw_legend:
            legend = self._ax.get_legend()
            if legend is not None:
                legend.remove()


def plot_channel_slice_grid(
    rows,
    plotter_cfg,
    out_path,
    want_data,
    signal_scale,
    title,
):
    """One figure per base category: slices across, lepton channels down.

    `rows` is [(channel, [panel or None, ...]), ...], each list holding one entry per
    slice in slice order. A missing slice keeps its column slot as None instead of
    shifting the row left -- dnn0 has to sit above dnn0 across channels, since making
    the columns comparable is the whole reason for sharing a canvas.

    Each panel keeps its own legend and channel label: the signal entry carries that
    panel's yield, so it is per-panel information rather than boilerplate.
    """
    import matplotlib

    if matplotlib.get_backend().lower() not in ("agg", "pdf", "svg", "ps"):
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import mplhep as hep

    from FLAF.PlotKit.plotters.stacked import StackedPlotter
    from FLAF.PlotKit.rootcompat import tlatex_to_mpl as _to_mathtext

    plt.style.use(hep.style.CMS)
    n_rows = len(rows)
    n_cols = max(len(panels) for _, panels in rows)
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(7.0 * n_cols, 7.5 * n_rows), squeeze=False
    )

    lumi_text = plotter_cfg.text_box("lumi_text").get("text", "")
    drawn_cms_label = False
    for r, (channel, panels) in enumerate(rows):
        for c in range(n_cols):
            ax = axes[r][c]
            panel = panels[c] if c < len(panels) else None
            if panel is None:
                # HistRebinTask produced no shapes here (e.g. a channel whose only
                # backgrounds are excluded). Blank the cell rather than closing the gap.
                ax.set_axis_off()
                continue
            slice_label, hists, custom = panel
            backend = PanelBackend(
                ax, draw_cms_label=not drawn_cms_label, draw_legend=True
            )
            drawn_cms_label = True
            sp = StackedPlotter(plotter_cfg, backend._impl)
            spec = sp.build_spec("rebinned", hists, want_data, custom, signal_scale)
            # The lumi belongs to the whole canvas, and inside a single narrow panel it
            # collides with "CMS Simulation" -- draw it once at figure level instead.
            spec.lumi_text = ""
            backend.render(spec)
            # Right-aligned: the CMS label the first panel draws is left-aligned above
            # the axes, and the lumi that would otherwise sit on the right was just
            # blanked, so this is the only thing on that side.
            ax.set_title(slice_label, loc="right", fontsize=16)
            if c > 0:
                ax.set_ylabel("")

    fig.suptitle(title, y=1.0, fontsize=20)
    if lumi_text:
        fig.text(
            1.0, 1.0, _to_mathtext(lumi_text), ha="right", va="bottom", fontsize=18
        )
    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def sum_over_eras(files, hist_path):
    total = None
    for f in files:
        h = f.Get(hist_path)
        if not h:
            continue
        h.SetDirectory(0)
        if total is None:
            total = h.Clone()
            total.SetDirectory(0)
        else:
            total.Add(h)
    return total


def build_category_panel(files, channel, category, cfg, mass, styles, era_label_cfg):
    """(hists, custom) for one slice, or None if the slice has no shapes.

    A slice is missing whenever HistRebinTask skipped the whole base category for
    lack of discovery signal (e.g. boosted at MX=300)."""
    signal_keys = [
        applyParameters(pattern, {cfg["signal_param_name"]: mass})
        for pattern in cfg["signal_hist_name_patterns"]
    ]
    prefix = f"{channel}/{category}/"

    hists = {}
    for bkg_key, allowed_channels in cfg["background_entries"]:
        if allowed_channels and channel not in allowed_channels:
            continue
        h = sum_over_eras(files, prefix + bkg_key)
        if h is None:
            continue
        label, color = styles.get(bkg_key, (bkg_key, "kGray"))
        hists[bkg_key] = (h, label, color, "backgrounds")

    if not hists:
        return None

    for signal_key in signal_keys:
        sig = sum_over_eras(files, prefix + signal_key)
        if sig is not None and sig.Integral() > 0:
            label, color = styles.get(signal_key, (signal_key, FALLBACK_SIGNAL_COLOR))
            hists[signal_key] = (sig, label, color, "signals")

    # category is "SR/res2b_dnn0"; the region prefix is already in the label boxes
    custom = {
        "cat_text": category.split("/")[-1],
        "ch_text": era_label_cfg.get("channel_text", {}).get(channel, channel),
        "customreg_text": era_label_cfg.get("customregion_text", {}).get(
            category.split("/")[0], ""
        ),
        # PlotKit already draws "Simulation" itself when want_data is False, so the
        # scope must not repeat it -- HistPlotter passes "CMS simulation" here and
        # ends up rendering "CMS Simulation simulation". Leaving the scope empty
        # defers to scope_text in the era's plot config.
        "datasim_text": "CMS ",
        "scope_text": "",
    }
    return hists, custom


def slice_title(files, channel, category):
    """Panel title: the sliced category name plus the selection it stands for.

    The selection is the slice directory's own title, set by hist_rebin_2d.py. Every era
    holds the same edges -- they are discovered once from the summed discovery eras -- so
    the first file carrying the slice wins. Input that was never sliced (a configuration
    with no `binning:` block) has no title of its own and simply carries no label.
    """
    name = category.split("/")[-1]
    for f in files:
        slice_dir = f.Get(f"{channel}/{category}")
        if not slice_dir:
            continue
        # ROOT defaults a directory's title to its own name, so only a title that
        # differs from it is a selection hist_rebin_2d.py actually wrote.
        label = slice_dir.GetTitle()
        return f"{name}   {label}" if label and label != name else name
    return name


def group_categories_by_base(categories, naming):
    """{"SR/res2b": ["SR/res2b_dnn0", ...]} in slice order -- the slices of one base
    category are what share a canvas."""
    groups = {}
    for cat in categories:
        groups.setdefault(naming.base(cat), []).append(cat)
    return groups


def run(
    input_dir,
    output_dir,
    config_path,
    eras,
    era_label,
    ana_path,
    version,
    signal_scale,
    log_y,
    masses=None,
    n_slices=None,
    category_pattern=None,
):
    from FLAF.PlotKit.config import PlotConfig

    cfg = load_config(config_path, n_slices, category_pattern)
    model = cfg["model"]
    mass_values = masses if masses else cfg["mass_values"]

    page_cfg = os.path.join(ana_path, "config", "plot", "cms_stacked.yaml")
    page_cfg_custom = os.path.join(ana_path, "config", "plot", f"{era_label}.yaml")
    if not os.path.exists(page_cfg_custom):
        print(f"Warning: {page_cfg_custom} not found, falling back to {eras[0]}")
        page_cfg_custom = os.path.join(ana_path, "config", "plot", f"{eras[0]}.yaml")
    with open(page_cfg_custom, "r") as f:
        era_label_cfg = yaml.safe_load(f)

    process_names = [name for name, _ in cfg["background_entries"]]
    process_names += [
        applyParameters(pattern, {cfg["signal_param_name"]: m})
        for pattern in cfg["signal_hist_name_patterns"]
        for m in mass_values
    ]
    styles = get_process_styles(ana_path, eras[0], version, process_names)

    # Every mass point's 2D variable projects onto the same y axis, so one lookup does
    # for all of them; getInputFileName() lays the tree out as "<era>/<variable>/...".
    first_rel = model.getInputFileName(
        eras[0], {cfg["signal_param_name"]: mass_values[0]}
    )
    variable = os.path.basename(os.path.dirname(first_rel))
    hist_cfg = build_hist_cfg(
        x_title=get_rebinned_axis_title(ana_path, eras[0], version, variable),
        y_title="Events",
        log_y=log_y,
    )
    plotter_cfg = PlotConfig(page_cfg, page_cfg_custom, hist_cfg)
    groups = group_categories_by_base(cfg["categories"], cfg["naming"])

    n_made, n_skipped = 0, 0
    for mass in mass_values:
        rel = model.getInputFileName(eras[0], {cfg["signal_param_name"]: mass})
        # getInputFileName() prefixes the era; input_dir holds one such tree per era
        var_rel = rel[len(eras[0]) + 1 :] if rel.startswith(eras[0] + "/") else rel

        files = []
        for era in eras:
            path = os.path.join(input_dir, era, var_rel)
            f = ROOT.TFile.Open(path, "READ")
            if f is None or f.IsZombie():
                raise RuntimeError(f"Cannot open rebinned file {path}")
            files.append(f)

        for base_cat, slice_cats in groups.items():
            # All lepton channels of a base category share one canvas, so the
            # slices can be compared across channels at a glance instead of by
            # flipping between files.
            rows = []
            for channel in cfg["channels"]:
                panels = []
                for category in slice_cats:
                    built = build_category_panel(
                        files, channel, category, cfg, mass, styles, era_label_cfg
                    )
                    if built is None:
                        panels.append(None)
                        continue
                    hists, custom = built
                    panels.append(
                        (
                            slice_title(files, channel, category),
                            hists,
                            custom,
                        )
                    )
                rows.append((channel, panels))

            if not any(panel for _, panels in rows for panel in panels):
                # HistRebinTask skipped this base category outright (e.g. boosted
                # at low MX): no slices exist in any channel, so there is nothing to draw.
                n_skipped += 1
                continue

            out_path = os.path.join(output_dir, f"m{mass}", f"{base_cat}.pdf")
            plot_channel_slice_grid(
                rows,
                plotter_cfg,
                out_path,
                want_data=False,
                signal_scale=signal_scale,
                title=f"{base_cat}   $M_X$ = {mass} GeV",
            )
            n_made += 1

        for f in files:
            f.Close()

    print(
        f"plot_rebinned: {n_made} canvases written to {output_dir} ({n_skipped} empty categories skipped)"
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Stacked plots of the merged, rebinned shapes entering the datacards."
    )
    parser.add_argument(
        "--input",
        required=True,
        type=str,
        help="base directory holding <era>/<var>/<var>.root for each era in --eras",
    )
    parser.add_argument(
        "--output", required=True, type=str, help="output directory for the pdfs"
    )
    parser.add_argument(
        "--config", required=True, type=str, help="datacard configuration yaml"
    )
    parser.add_argument(
        "--eras",
        required=True,
        type=str,
        help="comma-separated eras to sum before plotting",
    )
    parser.add_argument(
        "--era-label",
        required=True,
        type=str,
        help="name used for config/plot/<label>.yaml (lumi/channel/region text)",
    )
    parser.add_argument("--ana-path", required=True, type=str)
    parser.add_argument("--version", required=True, type=str)
    parser.add_argument(
        "--signal-scale",
        required=False,
        type=str,
        default="bkg",
        help="'bkg' to normalise signal to the summed background, or a fixed factor",
    )
    parser.add_argument("--log-y", action="store_true", help="log scale on y")
    parser.add_argument(
        "--masses",
        required=False,
        type=str,
        default=None,
        help="comma-separated subset of masses (default: all in the config)",
    )
    parser.add_argument(
        "--n-slices",
        required=False,
        type=int,
        default=None,
        help="slices each base category was cut into by HistRebinTask; defaults to "
        "the config's binning block. Pass the value that task actually used, so the "
        "panels match the input files",
    )
    parser.add_argument(
        "--category-pattern",
        required=False,
        type=str,
        default=None,
        help="pattern HistRebinTask named those slices with; defaults to the config's "
        "binning block. Pass the value that task actually used, so the panels carry the "
        "names the input files do",
    )
    args = parser.parse_args()

    try:
        signal_scale = float(args.signal_scale)
    except ValueError:
        signal_scale = args.signal_scale

    run(
        args.input,
        args.output,
        args.config,
        args.eras.split(","),
        args.era_label,
        args.ana_path,
        args.version,
        signal_scale,
        args.log_y,
        masses=args.masses.split(",") if args.masses else None,
        n_slices=args.n_slices,
        category_pattern=args.category_pattern,
    )
