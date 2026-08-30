import contextlib
import law
import luigi
import os

from string import Template

from FLAF.RunKit.run_tools import ps_call
from FLAF.run_tools.law_customizations import HTCondorWorkflow, copy_param

from StatInference.common.tools import CategoryNaming, importROOT

from .PreprocessShapesTask import PreprocessShapesTask
from .StatInferenceTask import StatInferenceTask


class CreateDatacardsTask(StatInferenceTask, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 2.0)
    n_cpus = copy_param(HTCondorWorkflow.n_cpus, 1)

    # If set, build datacards for this meta-era instead of self.period. self.period must
    # then be one of its real sub-eras -- it is only used to construct a valid FLAF Setup
    # (Setup.getGlobal requires a real, known period; a meta-era name isn't one).
    meta_era = luigi.Parameter(default="")

    # Stacked plots of the shapes going into the cards. This task is where they belong:
    # it is the only one holding every sub-era's histograms at once, which is what the
    # merged shape (and hence each datacard bin) is built from.
    make_plots = luigi.BoolParameter(default=True)

    @property
    def datacard_era(self):
        """The era datacards are actually built for: self.meta_era if set, else self.period."""
        return self.meta_era or self.period

    def get_sub_periods(self):
        """Real periods whose histograms feed into datacard_era: its constituent
        sub-eras for a meta-era, otherwise just [self.period]."""
        return self.get_era_groups().get(self.datacard_era, [self.period])

    def input_hist_reqs(self):
        """{key: task} for the shapes the cards are built from.

        The configuration's `preprocess:` step if it declares one, otherwise the merged
        histograms directly. Either way what arrives is "<era>/<variable>/<variable>.root",
        so nothing below here knows which it got.
        """
        if self.preprocess_config():
            return {
                "PreprocessShapes": PreprocessShapesTask.req(
                    self, meta_era=self.meta_era, branches=()
                )
            }
        return {
            f"MergedHists_{era}_{variable}": req
            for (era, variable), req in self.merged_hist_reqs(
                self.get_sub_periods()
            ).items()
        }

    def workflow_requires(self):
        return self.input_hist_reqs()

    def requires(self):
        return list(self.input_hist_reqs().values())

    def create_branch_map(self):
        return {0: None}

    def output(self):
        # fs_default. Note that combine cannot read these directly:
        # ResonantLimitsTask mirrors them back to datacards_dir() before handing them to
        # dhi -- see ResonantLimitsTask.stage_datacards.
        return self.output_dir_target(self.version, "Datacards", self.datacard_era)

    def run(self):
        statInf_entry = self.global_params["StatInference"]
        config = self.datacard_config_path()
        # ${ERA} in hist_bins names the era the cards are for, the same way the model's
        # input_file_pattern does. Each era has its own binning -- derived from its own
        # statistics, or its members' summed -- so one configuration serves them all.
        hist_bins_rel = statInf_entry.get("hist_bins")
        hist_bins = (
            os.path.join(
                self.ana_path(),
                Template(hist_bins_rel).safe_substitute(ERA=self.datacard_era),
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
                reqs = self.input_hist_reqs()["PreprocessShapes"]
                base_dir_local = stack.enter_context(reqs.output().localize("r"))
            else:
                targets = {
                    key: req.output()
                    for key, req in self.merged_hist_reqs(
                        self.get_sub_periods()
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
                self.datacard_era,
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

    def plot_variable(self, variable):
        """The variable whose axis the shapes are binned along, for histograms.yaml.

        For input a 2D->1D rebinning produced, that is the y variable of the 2D entry --
        histograms.yaml records both as ``var_list: [x, y]``, so the axis metadata the
        plotter needs is already described and does not have to be restated here. For
        input that was always 1D, the variable is its own answer.
        """
        try:
            import FLAF.Common.Setup as Setup

            hists = Setup.Setup(self.ana_path(), self.period, self.version).hists
            var_list = hists[variable].get("var_list")
            if var_list and len(var_list) > 1:
                return var_list[1]
        except Exception as e:
            print(f"Warning: no var_list for {variable} ({e}); plotting it as itself")
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
    def _present_keys(shape_file, cfg):
        """The (channel, region, category) triples the summed shape file actually holds."""
        ROOT = importROOT()
        present = set()
        f = ROOT.TFile.Open(shape_file)
        if not f or f.IsZombie():
            return present
        try:
            for category in cfg["categories"]:
                region, _, cat = category.rpartition("/")
                for channel in cfg["channels"]:
                    path = "/".join(x for x in (channel, region, cat) if x)
                    if f.Get(path):
                        present.add((channel, region, cat))
        finally:
            f.Close()
        return present

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
        plots_dir = os.path.join(output_dir, "plots")
        os.makedirs(plots_dir, exist_ok=True)
        plotter = os.path.join(self.ana_path(), "FLAF", "Analysis", "HistPlotter.py")

        # The datacard bin is the sum over the sub-eras (getCombinedShape does the same at
        # card-build time), so the sub-era files are hadd'ed into the one file the plotter
        # reads. Plotting them separately would show something the fit never sees.
        for variable in self.get_required_variables():
            inputs = [
                p
                for era in self.get_sub_periods()
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
            present = self._present_keys(summed, cfg)
            keys, outputs = [], []
            for category in cfg["categories"]:
                # "SR/res2b_dnn0" -> region "SR", category "res2b_dnn0": HistPlotter
                # navigates channel -> region -> category.
                region, _, cat = category.rpartition("/")
                for channel in cfg["channels"]:
                    if (channel, region, cat) not in present:
                        continue
                    keys.append(f"{channel}:{cat}:{region}")
                    outputs.append(
                        os.path.join(
                            plots_dir, f"{variable}_{channel}_{region}_{cat}.pdf"
                        )
                    )
            if not keys:
                print(f"WARNING: no shapes to plot for {variable}")
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
                "--globalConfig",
                os.path.join(
                    self.ana_path(),
                    self.global_params["analysis_config_area"],
                    "global.yaml",
                ),
                # The surviving axis of the rebinned shapes, and already at its final
                # binning -- so no --rebin, which would coarsen it back to the histograms.yaml
                # grid and undo the whole point of the rebinning.
                "--var",
                self.plot_variable(variable),
                # The plotted shapes are the sum over every sub-era, so the label has to
                # name the combination: HistPlotter reads config/plot/<year>.yaml for the
                # luminosity, and a sub-era's file states that sub-era's luminosity alone.
                # --period stays a real era -- it builds the Setup, which does not know
                # meta-era names.
                "--year",
                self.datacard_era,
                "--ana_path",
                self.ana_path(),
                "--period",
                self.period,
                "--LAWrunVersion",
                self.version,
                "--wantSignals",
                # The slices span ~5 decades in yield (the low-significance slice holds
                # most of the background), so a linear axis hides everything but slice 0.
                "--wantLogScale",
                "y",
            ]
            try:
                ps_call(cmd, verbose=1)
                self._stitch_grids(cfg, plots_dir, variable, dict(zip(keys, outputs)))
            except Exception as e:
                # The datacards are the product that matters; a plotting failure should be
                # loud but must not leave the task looking broken with valid cards on disk.
                print(
                    f"WARNING: shape plotting failed for {variable}, "
                    f"datacards are unaffected: {e}"
                )
            finally:
                if os.path.exists(summed):
                    os.remove(summed)
