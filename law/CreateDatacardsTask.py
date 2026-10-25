import contextlib
import law
import luigi
import os

from string import Template

from FLAF.RunKit.run_tools import ps_call
from FLAF.run_tools.law_customizations import (
    EraDependentValueError,
    HTCondorWorkflow,
    copy_param,
)

from StatInference.common.tools import CategoryNaming, InputCategories, importROOT

from .PreprocessShapesTask import PreprocessShapesTask
from .StatInferenceTask import StatInferenceTask


class CreateDatacardsTask(StatInferenceTask, HTCondorWorkflow, law.LocalWorkflow):
    """The datacards of the configuration, one branch per card (`eras:`, in its order).

    Runs dc_make/create_datacards.py on whatever `<era>/<variable>/<variable>.root` tree
    it is given -- the configuration's `preprocess:` step if it declares one, otherwise
    the merged histograms directly -- and does no binning of its own. The card of a group
    era holds every sub-era's histograms at once, which is what the summed shape behind
    each datacard bin is built from, and why the stacked shape plots are drawn here.
    """

    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 2.0)
    n_cpus = copy_param(HTCondorWorkflow.n_cpus, 1)

    # Stacked plots of the shapes going into the cards. This task is where they belong:
    # it is the only one holding every sub-era's histograms at once, which is what the
    # merged shape (and hence each datacard bin) is built from.
    make_plots = luigi.BoolParameter(default=True)

    def input_hist_reqs(self, card_eras):
        """{key: task} for the shapes the cards of `card_eras` are built from.

        The configuration's `preprocess:` step if it declares one -- its workflow over the
        same cards (law hands the --branches selection on), or from a branch its branch of
        the same card -- otherwise the merged histograms directly. Either way what arrives
        is "<era>/<variable>/<variable>.root", so nothing below here knows which it got.
        """
        if self.preprocess_config():
            return {"PreprocessShapes": PreprocessShapesTask.req(self)}
        return {
            f"MergedHists_{era}_{variable}": req
            for (era, variable), req in self.merged_hist_reqs(
                self.sub_eras_of(card_eras)
            ).items()
        }

    def stage_workflow_requires(self):
        # Merged with the base class's rather than replacing them: FLAF's HTCondorWorkflow
        # puts the software bundles a submitted job unpacks in there, and returning only
        # our own inputs drops them, so the jobs start without the bundle they need.
        reqs = super().stage_workflow_requires()
        # the selected cards only (the map honours --branches)
        reqs.update(self.input_hist_reqs(self.get_branch_map().values()))
        return reqs

    def requires(self):
        return list(self.input_hist_reqs([self.branch_data]).values())

    def create_branch_map(self):
        return dict(enumerate(self.card_eras()))

    def output(self):
        # fs_default. Note that combine cannot read these directly:
        # ResonantLimitsTask mirrors them back to datacards_dir() before handing them to
        # dhi -- see ResonantLimitsTask.stage_datacards.
        return self.output_dir_target(self.version, "Datacards", self.branch_data)

    def run(self):
        card_era = self.branch_data
        statInf_entry = self.statinference_settings()
        config = self.datacard_config_path()
        # ${ERA} in hist_bins names the era the cards are for, the same way the model's
        # input_file_pattern does. Each era has its own binning -- derived from its own
        # statistics, or its members' summed -- so one configuration serves them all.
        hist_bins_rel = statInf_entry.get("hist_bins")
        hist_bins = (
            os.path.join(
                self.ana_path(),
                Template(hist_bins_rel).safe_substitute(ERA=card_era),
            )
            if hist_bins_rel
            else None
        )
        param_values = statInf_entry.get("param_values", [])
        create_datacards_py = os.path.join(
            self.ana_path(), "StatInference", "dc_make", "create_datacards.py"
        )
        with contextlib.ExitStack() as stack:
            if self.preprocess_config():
                # PreprocessShapesTask already wrote the "<era>/<variable>/<variable>.root"
                # layout, so its output directory is the base input_file_pattern resolves
                # against -- nothing to stage.
                reqs = self.input_hist_reqs([card_era])["PreprocessShapes"]
                base_dir_local = stack.enter_context(reqs.output().localize("r"))
            else:
                targets = {
                    key: req.output()
                    for key, req in self.merged_hist_reqs(
                        self.get_sub_eras(card_era)
                    ).items()
                }
                self.check_inputs(targets)
                base_dir_local = stack.enter_context(self.stage_inputs(targets))
            local_output = stack.enter_context(self.output().localize("w"))
            cmd = [
                "python3",
                create_datacards_py,
                "--input",
                base_dir_local.abspath,
                "--output",
                local_output.abspath,
                "--config",
                config,
                "--eras",
                card_era,
            ]
            if hist_bins:
                cmd += ["--hist-bins", hist_bins]
            if len(param_values) > 0:
                param_values_str = ",".join(str(v) for v in param_values)
                cmd += ["--param_values", param_values_str]
            ps_call(cmd, env=self.cmssw_env, verbose=1)

            if self.make_plots:
                self.plot_rebinned_shapes(
                    base_dir_local.abspath, config, local_output.abspath
                )

    @staticmethod
    def _surviving_axis(base_dir):
        """Which axis of the 2D input the rebinned shapes are binned along: 0 = x, 1 = y.

        Read from binning.json: a slice with a "y_range" (a window on y) is binned along x,
        anything else along y. Without a record, y.
        """
        import json

        path = os.path.join(base_dir, "binning.json")
        try:
            with open(path) as f:
                record = json.load(f)
        except (OSError, ValueError):
            return 1
        node = record.get("binning", {})
        # era -> mass -> channel -> category -> {"slices": [...]}
        while isinstance(node, dict) and "slices" not in node:
            if not node:
                return 1
            node = next(iter(node.values()))
        slices = node.get("slices") if isinstance(node, dict) else None
        return 0 if slices and "y_range" in slices[0] else 1

    def plot_variable(self, variable, sub_eras, axis=1):
        """The variable whose axis the shapes are binned along, for histograms.yaml.

        For 2D->1D rebinned input, the entry of the 2D variable's ``var_list: [x, y]``
        picked by `axis` (from _surviving_axis()); for 1D input, the variable itself. The
        plotted shapes are the sum over `sub_eras`, so their histograms.yaml must agree.
        """
        # None: the era's histograms.yaml has no entry for the variable
        var_lists = {}
        for era in sub_eras:
            hists = self.setup_for(era).hists
            try:
                var_lists[era] = hists[variable].get("var_list")
            except KeyError:
                var_lists[era] = None
        var_list = var_lists[sub_eras[0]]
        if any(other != var_list for other in var_lists.values()):
            raise EraDependentValueError(
                f"{self.get_task_family()}: the var_list of {variable} in histograms.yaml "
                f"differs between the eras summed in one plot ({var_lists})"
            )
        if variable not in self.setup_for(sub_eras[0]).hists:
            print(
                f"Warning: no histograms.yaml entry for {variable}; plotting it as itself"
            )
        elif var_list and len(var_list) > 1:
            return var_list[axis]
        return variable

    @staticmethod
    def _stitch_grid(panels, out_path):
        """Lay rendered panels out on one page: a row per channel, a column per slice.

        The slices of a base category are one physical selection cut into pieces, and the
        fit sees them together -- a slice that looks reasonable in eMu and pathological in
        eE is obvious side by side and invisible in separate files. Only page geometry
        happens here; every panel is drawn by HistPlotter.py, so there is still one
        plotting implementation.

        `panels` is [[path or None, ...] per column] per row; a None leaves its cell blank,
        which is what keeps column N under column N when a slice was skipped.
        """
        from pypdf import PageObject, PdfReader, PdfWriter, Transformation

        pages = {
            p: PdfReader(p).pages[0] for row in panels for p in row if p is not None
        }
        if not pages:
            return False
        w = max(float(p.mediabox.width) for p in pages.values())
        h = max(float(p.mediabox.height) for p in pages.values())
        n_rows, n_cols = len(panels), max(len(r) for r in panels)

        sheet = PageObject.create_blank_page(width=w * n_cols, height=h * n_rows)
        for r, row in enumerate(panels):
            for c, path in enumerate(row):
                if path is None:
                    continue
                # PDF origin is bottom-left, so the first row goes at the top.
                sheet.merge_transformed_page(
                    pages[path],
                    Transformation().translate(c * w, (n_rows - 1 - r) * h),
                )
        writer = PdfWriter()
        writer.add_page(sheet)
        with open(out_path, "wb") as f:
            writer.write(f)
        return True

    @staticmethod
    def _present_keys(shape_file, cfg, points):
        """({(channel, region, category): input directory} for the datacard bins the summed
        shape file holds, [every path looked up], or None when the file cannot be read).

        Each category is looked up where the maker reads it for these model points: its
        `input_categories` directory, else its own name. A category that the points read
        from different directories is not plotted, since one plot cannot show them.
        """
        ROOT = importROOT()
        input_categories = InputCategories.fromConfig(cfg)
        present, looked_up = {}, []
        f = ROOT.TFile.Open(shape_file)
        if not f or f.IsZombie():
            return present, None
        try:
            for category in cfg["categories"]:
                region, _, cat = category.rpartition("/")
                directories = {input_categories.path(category, p) for p in points}
                if len(directories) != 1:
                    print(
                        f"WARNING: {category} is read from {', '.join(sorted(directories))} "
                        "for the model points of one variable; it is not plotted"
                    )
                    continue
                (directory,) = directories
                for channel in cfg["channels"]:
                    path = "/".join(x for x in (channel, directory) if x)
                    looked_up.append(path)
                    if f.Get(path):
                        present[(channel, region, cat)] = directory
        finally:
            f.Close()
        return present, looked_up

    def _stitch_grids(self, cfg, plots_dir, variable, panel_of_key):
        """One grid per base category: its slices across, the channels down."""
        naming = CategoryNaming.fromConfig(cfg)
        bases = {}
        for category in cfg["categories"]:
            region, _, cat = category.rpartition("/")
            base, slice_idx = naming.split(cat)
            bases.setdefault((region, base), {})[slice_idx] = cat
        for (region, base), slices in bases.items():
            # None for an unsliced category, which is a single-column grid of channels.
            order = sorted(slices, key=lambda i: (i is None, i))
            panels = [
                [
                    (lambda p: p if p and os.path.exists(p) else None)(
                        panel_of_key.get(f"{channel}:{slices[i]}:{region}")
                    )
                    for i in order
                ]
                for channel in cfg["channels"]
            ]
            out = os.path.join(
                plots_dir, f"{variable}_{region.replace('/', '_')}_{base}_grid.pdf"
            )
            if self._stitch_grid(panels, out):
                print(f"Wrote grid {out}")

    def plot_rebinned_shapes(self, base_dir, config, output_dir):
        """Stacked plots of the shapes the cards are built from, one per datacard bin.

        Runs FLAF's own HistPlotter.py -- the script HistPlotTask uses -- rather than a
        plotter of our own, so these come out in the same style as every other plot in the
        analysis and there is one implementation to maintain. It is agnostic about which
        categories exist: it plots whatever `channel:category:region` keys it is handed,
        so the sliced names need nothing on the FLAF side.

        Runs in the default env rather than cmssw_env: the plotter needs PlotKit
        (matplotlib/mplhep), same as HistPlotTask.
        """
        import glob

        cfg = self.get_config_data()
        card_era = self.branch_data
        sub_eras = self.get_sub_eras(card_era)
        plots_dir = os.path.join(output_dir, "plots")
        os.makedirs(plots_dir, exist_ok=True)
        plotter = os.path.join(self.ana_path(), "FLAF", "Analysis", "HistPlotter.py")

        # The datacard bin is the sum over the sub-eras (getCombinedShape does the same at
        # card-build time), so the sub-era files are hadd'ed into the one file the plotter
        # reads. Plotting them separately would show something the fit never sees.
        variable_points = self.get_variable_points()
        for variable in self.get_required_variables():
            inputs = [
                p
                for era in sub_eras
                for p in glob.glob(
                    os.path.join(base_dir, era, variable, f"{variable}.root")
                )
            ]
            if not inputs:
                continue
            summed = os.path.join(plots_dir, f"_summed_{variable}.root")
            try:
                ps_call(["hadd", "-f", summed] + inputs, verbose=1)
            except Exception as e:
                print(
                    f"WARNING: could not sum the shapes for {variable}, "
                    f"datacards are unaffected: {e}"
                )
                continue

            # Only the categories this mass point actually has. The configuration lists
            # every datacard bin, but the preprocessing gates drop a category where the
            # signal is too small -- boosted at a low mass, say -- and it is then absent
            # from the shapes. HistPlotter exits non-zero on the first key it cannot find,
            # so passing the full list costs every panel and grid that would have come
            # after the gap, not just the missing one.
            present, looked_up = self._present_keys(
                summed, cfg, variable_points[variable]
            )
            keys, outputs, panel_of_key = [], [], {}
            for category in cfg["categories"]:
                # "SR/res2b_dnn0" -> region "SR", category "res2b_dnn0": HistPlotter
                # navigates channel -> region -> category.
                region, _, cat = category.rpartition("/")
                for channel in cfg["channels"]:
                    directory = present.get((channel, region, cat))
                    if directory is None:
                        continue
                    # The key names the directory the shapes are read from; the file and
                    # the grid keep the datacard bin's name.
                    in_region, _, in_cat = directory.rpartition("/")
                    keys.append(f"{channel}:{in_cat}:{in_region}")
                    output = os.path.join(
                        plots_dir, f"{variable}_{channel}_{region}_{cat}.pdf"
                    )
                    outputs.append(output)
                    panel_of_key[f"{channel}:{cat}:{region}"] = output
            if not keys:
                if looked_up is None:
                    reason = f"{summed} cannot be read"
                elif looked_up:
                    reason = f"none of {', '.join(looked_up)} is in the summed shapes"
                else:
                    reason = "no category is read from one directory"
                print(f"WARNING: no shapes to plot for {variable}: {reason}")
                if os.path.exists(summed):
                    os.remove(summed)
                continue
            cmd = [
                "python3",
                plotter,
                "--inFile",
                summed,
                "--all_outFiles",
                ",".join(outputs),
                "--all_keys",
                ",".join(keys),
                # No config file: HistPlotter takes every setting from the Setup it builds
                # from --ana_path, --era and --user-custom, as the task's own Setup does.
                # The surviving axis of the rebinned shapes, and already at its final
                # binning -- so no --rebin, which would coarsen it back to the histograms.yaml
                # grid and undo the whole point of the rebinning.
                "--var",
                self.plot_variable(variable, sub_eras, self._surviving_axis(base_dir)),
                # The plotted shapes are the sum over every sub-era, so the label has to
                # name the combination: HistPlotter reads config/plot/<year>.yaml for the
                # luminosity, and a sub-era's file states that sub-era's luminosity alone.
                # --era stays a real era -- it builds the Setup, which does not know
                # meta-era names.
                "--year",
                card_era,
                "--ana_path",
                self.ana_path(),
                "--era",
                sub_eras[0],
                "--LAWrunVersion",
                self.version,
                "--wantSignals",
                # The slices span ~5 decades in yield (the low-significance slice holds
                # most of the background), so a linear axis hides everything but slice 0.
                "--wantLogScale",
                "y",
            ]
            if self.user_custom:
                cmd += ["--user-custom", self.user_custom]
            try:
                ps_call(cmd, verbose=1)
                self._stitch_grids(cfg, plots_dir, variable, panel_of_key)
            except Exception as e:
                # A failed plot fails the task: a warning here hid a broken HistPlotter call
                # behind a green CI (FLAF#365). The cards themselves are fine, so say how to
                # make them without the plots.
                raise RuntimeError(
                    f"shape plotting failed for {variable}: {e}. The datacards are "
                    "unaffected; run with --make-plots False to produce them without plots."
                ) from e
            finally:
                if os.path.exists(summed):
                    os.remove(summed)
