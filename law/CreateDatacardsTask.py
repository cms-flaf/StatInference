import contextlib
import law
import luigi
import os

from FLAF.RunKit.run_tools import ps_call
from FLAF.run_tools.law_customizations import HTCondorWorkflow, copy_param

from .RebinningParamsTask import RebinningParamsTask
from .HistRebinTask import HistRebinTask


class CreateDatacardsTask(RebinningParamsTask, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 2.0)
    n_cpus = copy_param(HTCondorWorkflow.n_cpus, 1)

    # If set, build datacards for this meta-era instead of self.period. self.period must
    # then be one of its real sub-eras -- it is only used to construct a valid FLAF Setup
    # (Setup.getGlobal requires a real, known period; a meta-era name isn't one).
    meta_era = luigi.Parameter(default="")

    # Stacked plots of the rebinned shapes going into the cards. This task is where they
    # belong: it is the only one holding every sub-era's HistRebinTask output at once,
    # which is what the merged shape (and hence each datacard bin) is built from.
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

        With rebinning, that is one HistRebinTask per sub-period; without it (a
        configuration whose input is already 1D and binned) the merged histograms are read
        straight from the Hists_merged tree, since there is no step in between.
        """
        if not self.rebinning_enabled():
            return {
                f"MergedHists_{era}_{variable}": req
                for (era, variable), req in self.merged_hist_reqs(
                    self.get_sub_periods()
                ).items()
            }
        return {
            f"HistRebin_{e}": HistRebinTask.req(
                self, period=e, meta_era=self.meta_era, branches=()
            )
            for e in self.get_sub_periods()
        }

    def workflow_requires(self):
        return self.input_hist_reqs()

    def requires(self):
        return list(self.input_hist_reqs().values())

    def create_branch_map(self):
        return {0: None}

    def output(self):
        # fs_default, like HistRebinTask. Note that combine cannot read these directly:
        # ResonantLimitsTask mirrors them back to datacards_dir() before handing them to
        # dhi -- see ResonantLimitsTask.stage_datacards.
        return self.output_dir_target(
            self.version, "Datacards", *self.binning_parts(), self.datacard_era
        )

    def run(self):
        statInf_entry = self.global_params["StatInference"]
        config = self.datacard_config_path()
        hist_bins_rel = statInf_entry.get("hist_bins")
        hist_bins = (
            os.path.join(self.ana_path(), hist_bins_rel) if hist_bins_rel else None
        )
        param_values = statInf_entry.get("param_values", [])
        create_datacards_py = os.path.join(
            self.ana_path(), "StatInference", "dc_make", "create_datacards.py"
        )
        with contextlib.ExitStack() as stack:
            if self.rebinning_enabled():
                # HistRebinTask writes <era>/<variable>/<variable>.root under its own
                # per-period directory, so its parent is already the base the config's
                # input_file_pattern resolves against.
                base_dir_local = stack.enter_context(
                    self.input()[0].parent.localize("r")
                )
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
            if self.rebinning_enabled():
                # The datacard bins are the sliced category names; the resolved values
                # here are what HistRebinTask actually wrote, and must win over the
                # config's own so the two cannot drift apart under a CLI override.
                # Left unset otherwise, so the categories stay as the config lists them.
                binning = self.binning_params()
                cmd += ["--n-slices", str(binning["n_slices"])]
                cmd += ["--category-pattern", binning["category_pattern"]]
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

    def plot_rebinned_shapes(self, base_dir, config, output_dir):
        """Stacked plots of the sub-era-summed rebinned shapes, one per datacard bin.

        Runs in the default env rather than cmssw_env: the plotter needs FLAF's PlotKit
        (matplotlib/mplhep), same as HistPlotTask.
        """
        plot_py = os.path.join(
            self.ana_path(), "StatInference", "dc_make", "plot_rebinned.py"
        )
        cmd = [
            "python3",
            plot_py,
            "--input",
            base_dir,
            "--output",
            os.path.join(output_dir, "plots"),
            "--config",
            config,
            "--eras",
            ",".join(self.get_sub_periods()),
            "--era-label",
            self.datacard_era,
            "--ana-path",
            self.ana_path(),
            "--version",
            self.version,
            "--signal-scale",
            str(self.global_params.get("signal_plot_scale", "bkg")),
            # The slices span ~5 decades in yield (the low-significance slice holds
            # most of the background), so a linear axis hides everything but slice 0.
            "--log-y",
        ]
        if self.rebinning_enabled():
            # Same reason as create_datacards.py above: the panels are the sliced
            # categories, so they must follow the values HistRebinTask actually used.
            # Left unset otherwise, so the panels are the config's own categories.
            binning = self.binning_params()
            cmd += ["--n-slices", str(binning["n_slices"])]
            cmd += ["--category-pattern", binning["category_pattern"]]
        try:
            ps_call(cmd, verbose=1)
        except Exception as e:
            # The datacards are the product that matters; a plotting failure should be
            # loud but must not leave the task looking broken with valid cards on disk.
            print(
                f"WARNING: rebinned-shape plotting failed, datacards are unaffected: {e}"
            )
