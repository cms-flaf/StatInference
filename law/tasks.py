import contextlib
import glob
import hashlib
import json
import law
import luigi
import math
import os
import re
import shutil
import subprocess
import yaml

from string import Template

from FLAF.RunKit.run_tools import ps_call
from FLAF.run_tools.law_customizations import (
    Task,
    HTCondorWorkflow,
    copy_param,
)
from FLAF.Analysis.tasks import HistPlotTask
from StatInference.common.tools import CategoryNaming
from dhi.config import campaign_labels, campaign_lumis
from dhi.tasks.resonant import (
    MergeResonantLimits,
    PlotMultipleResonantLimits,
    PlotResonantLimits,
)


class StatInferenceTask(Task):
    """Shared datacard-configuration access for the limit-setting chain.

    Everything downstream of the merged histograms -- rebinning, datacards, limits, limit
    plots -- is driven by the datacard configuration named in global.yaml's
    ``StatInference.config``, not by global.yaml's own variable lists. This base class is
    the single place that file is read, so the four tasks below cannot disagree about
    which eras, masses or variables the analysis consists of.
    """

    # The Hists_merged tree to read. Separate from `version`, which names what this chain
    # *writes*: re-binning or re-fitting under a new version must not require the input
    # histograms to be reproduced (or copied) under that name. Not significant -- it
    # identifies an input, not a product, and law's store paths are about products.
    hists_version = luigi.Parameter(
        default="",
        significant=False,
        description="version of the Hists_merged tree to read; defaults to --version",
    )

    @property
    def input_hists_version(self):
        return self.hists_version or self.version

    def output_dir_target(self, *path):
        """remote_dir_target() that also works when fs_default is a local directory.

        FLAF's remote_dir_target() handles a str fs and a remote fs but not a
        law.LocalFileSystem: it falls through to WLCGDirectoryTarget(), which raises
        "fs must be a RemoteFileSystem instance". Its sibling remote_target() (the file
        variant) already carries the branch mirrored here. The CI configures fs_default
        as a plain local path, so every directory output in this chain has to cope --
        without this, output() raises and luigi reports the task as having "error in
        complete() method".
        """
        fs = self.fs_default
        if isinstance(fs, law.LocalFileSystem):
            return law.LocalDirectoryTarget(os.path.join(*path), fs=fs)
        return self.remote_dir_target(*path)

    def datacard_config_path(self):
        return os.path.join(
            self.ana_path(), self.global_params["StatInference"]["config"]
        )

    def get_config_data(self):
        # Cached per instance: requires()/workflow_requires() are re-entered many times
        # during graph construction and each call would otherwise re-parse the yaml.
        if getattr(self, "_config_data", None) is None:
            with open(self.datacard_config_path(), "r") as f:
                self._config_data = yaml.safe_load(f)
        return self._config_data

    def get_era_groups(self):
        return self.get_config_data().get("era_groups", {})

    def get_top_level_eras(self):
        """Eras that carry their own datacards and limits: meta-eras plus any standalone
        real era. Real eras already covered by a meta-era are excluded so they are not
        double-counted alongside it."""
        data = self.get_config_data()
        eras = data.get("eras", [self.period])
        grouped = {e for sub in self.get_era_groups().values() for e in sub}
        return [e for e in eras if e not in grouped]

    def get_campaign(self, era):
        """dhi campaign key for an era, from the config's ``campaigns`` map.

        It cannot be derived from the era name. 'Run3_Early' happens to lowercase into
        dhi's 'run3_early', but 'Run3_2022' would give 'run3_2022', which is not a key of
        dhi.config.campaign_labels -- and dhi draws an unknown key verbatim rather than
        failing, so a guessed value ends up printed on the plot as the luminosity label.
        """
        return self.get_config_data().get("campaigns", {}).get(era)

    def get_masses(self):
        """Every parameter value the configuration declares, across all processes.

        Not restricted to signals: with ``param_dependent_bkg`` a background can be
        parameterised too, and its histograms live in the same per-mass input file.
        """
        masses = set()
        for proc in self.get_config_data().get("processes", []):
            if not isinstance(proc, dict):
                # A bare string entry names a background directly and carries no
                # settings -- the shape dc_make's load_config functions also accept.
                continue
            for value in proc.get("param_values") or []:
                masses.add(value)
        if not masses:
            raise RuntimeError(
                f"No process in {self.datacard_config_path()} declares param_values, so "
                "there is no mass point to build inputs for."
            )
        return sorted(masses)

    def get_required_variables(self):
        """The Hists_merged variables this chain reads: one per mass point (the 2D
        DNN-vs-HME histogram, or the 1D DNN score for a configuration that does not
        rebin), derived from the datacard config's input_file_pattern.

        This is the *complete* input list -- there is no intersection with global.yaml's
        histTuple_flavor variable list anywhere downstream, so a mass point present here
        is read whether or not the active flavour happens to mention it.
        """
        data = self.get_config_data()
        model = data["model"]
        pattern = model["input_file_pattern"]
        param_name = model["parameters"][0]

        variables = set()
        for mass in self.get_masses():
            rel = Template(pattern).safe_substitute({"ERA": "ERA", param_name: mass})
            # "<era>/<variable>/<variable>.root": the variable is the directory holding
            # the file, not the file itself.
            variables.add(os.path.basename(os.path.dirname(rel)))
        return sorted(variables)

    def merged_hist_reqs(self, eras):
        """{(era, variable): MergedHists} -- every merged input file for those eras.

        Shared by the two tasks that read Hists_merged: HistRebinTask (over its discovery
        eras) and, for a configuration with no rebinning step, CreateDatacardsTask (over
        its sub-periods).
        """
        return {
            (era, variable): MergedHists.req(self, period=era, variable=variable)
            for era in eras
            for variable in self.get_required_variables()
        }

    @contextlib.contextmanager
    def stage_inputs(self, targets):
        """Assemble the required merged histograms into the "<era>/<variable>/
        <variable>.root" layout that Model.getInputFileName resolves input_file_pattern
        against, and yield the directory holding it.

        Each input is localized individually so only the files this task actually reads
        cross the network -- the merged tree holds every variable of the active
        histTuple_flavor (~106 for 'default') across every era.
        """
        with contextlib.ExitStack() as stack:
            staging = law.LocalDirectoryTarget(is_tmp=True)
            staging.touch()
            stack.callback(lambda: staging.remove(silent=True))

            for (era, variable), target in targets.items():
                dest_dir = os.path.join(staging.abspath, era, variable)
                os.makedirs(dest_dir, exist_ok=True)
                local_inp = stack.enter_context(target.localize("r"))
                shutil.copy2(
                    local_inp.abspath, os.path.join(dest_dir, f"{variable}.root")
                )

            yield staging

    def check_inputs(self, targets):
        """Report every missing input at once, naming the version they were looked for
        under. law already refuses to run with an incomplete MergedHists dependency; this
        makes the condor log self-explanatory when the task is forced anyway, instead of
        failing on whichever file ROOT happened to open first."""
        missing = sorted(
            f"{era}/{var}" for (era, var), t in targets.items() if not t.exists()
        )
        if missing:
            raise RuntimeError(
                f"{len(missing)} of {len(targets)} merged histograms are missing under "
                f"hists_version='{self.input_hists_version}' "
                f"(<fs_HistTuple>/{self.input_hists_version}/Hists_merged/<era>/<var>/<var>.root): "
                + ", ".join(missing)
            )


class MergedHists(StatInferenceTask, law.ExternalTask):
    """One merged histogram file, addressed by path under --hists-version.

    External on purpose. The limit-setting chain consumes histograms it does not make, so
    a missing input is reported as a missing dependency rather than pulling the whole
    AnaTuple production graph (AnaTupleFileListTask -> HistFromNtupleProducerTask ->
    HistMergerTask) into the run and rebuilding its branch maps just to discover the file
    is already there.

    The path mirrors HistMergerTask.output() in FLAF/Analysis/tasks.py -- that is the
    contract between the two halves, and the only thing this chain needs from FLAF.
    """

    variable = luigi.Parameter(description="Hists_merged variable (directory) name")

    def output(self):
        return self.remote_target(
            self.input_hists_version,
            "Hists_merged",
            self.period,
            self.variable,
            f"{self.variable}.root",
            fs=self.fs_HistTuple,
        )


# Binning knobs whose values change the rebinned shapes. Kept in one place because they
# are propagated down the chain, read from the datacard configuration's `binning:` block,
# and hashed into the output path -- see RebinningParamsTask.binning_parts().
BINNING_PARAMS = (
    "n_dnn_slices",
    "category_pattern",
    "slice_var",
    "max_mass_bins",
    "min_dnn_bkg_sum",
    "min_mass_bkg_each",
    "min_signal",
    "min_bkg_neff",
    "min_bkg_frac",
    "min_mass_bkg_neff",
    "bkg_per_mass_bin",
    "significance_mode",
    "min_dnn_bkg_each",
    "min_dnn_bkg_neff",
)

# Fallbacks for a configuration that declares no `binning:` block. The values that matter
# belong to an analysis, not to this file: what a slice needs to be worth keeping depends
# on the sample sizes and the selection, and a test configuration wants them at zero. Each
# datacard configuration should state its own -- see the annotated block in
# config/Datacards/x_hh_bbww_DL_run3.yaml for what the production numbers are and why.
BINNING_DEFAULTS = {
    "n_dnn_slices": 4,
    # Neutral by default: the sliced axis is whatever the analysis put on the 2D x axis,
    # and an analysis that wants its own name for it says so in its `binning:` block.
    "category_pattern": CategoryNaming.default_pattern,
    "slice_var": "x",
    "max_mass_bins": 10,
    "min_dnn_bkg_sum": 1.0,
    "min_mass_bkg_each": 0.01,
    "min_signal": 0.5,
    "min_bkg_neff": 4.0,
    "min_bkg_frac": 0.05,
    "min_mass_bkg_neff": 4.0,
    "bkg_per_mass_bin": 5.0,
    "significance_mode": "asimov",
    "min_dnn_bkg_each": 0.01,
    "min_dnn_bkg_neff": 0.0,
}


class RebinningParamsTask(StatInferenceTask):
    """The binning parameters, declared once for the whole chain.

    They live here rather than on HistRebinTask alone for two reasons. First, `.req()`
    only copies parameters the requesting task also declares, so with them on HistRebinTask
    only there was no way to steer the binning from ResonantLimitsTask -- the top of the
    chain always got the defaults. Second, every task's output path has to agree on which
    binning it describes, which needs the values everywhere.
    """

    # All default to None, meaning "take it from the datacard configuration". Optional*
    # rather than plain Int/Float/Parameter is required, not cosmetic: luigi serialises a
    # plain FloatParameter's None to the string "None" and then raises on float("None")
    # whenever it round-trips the task through to_str_params()/from_str_params().
    n_dnn_slices = luigi.OptionalIntParameter(default=None)
    category_pattern = luigi.OptionalParameter(default=None)
    slice_var = luigi.OptionalParameter(default=None)
    max_mass_bins = luigi.OptionalIntParameter(default=None)
    min_dnn_bkg_sum = luigi.OptionalFloatParameter(default=None)
    min_mass_bkg_each = luigi.OptionalFloatParameter(default=None)
    min_signal = luigi.OptionalFloatParameter(default=None)
    min_bkg_neff = luigi.OptionalFloatParameter(default=None)
    min_bkg_frac = luigi.OptionalFloatParameter(default=None)
    min_mass_bkg_neff = luigi.OptionalFloatParameter(default=None)
    bkg_per_mass_bin = luigi.OptionalFloatParameter(default=None)
    significance_mode = luigi.OptionalParameter(default=None)
    min_dnn_bkg_each = luigi.OptionalFloatParameter(default=None)
    min_dnn_bkg_neff = luigi.OptionalFloatParameter(default=None)

    def rebinning_enabled(self):
        """Whether the configuration asks for the in-chain 2D->1D rebinning step.

        A configuration with no ``binning:`` block describes input that is already 1D and
        binned: HistRebinTask has nothing to do with it, and CreateDatacardsTask reads the
        merged histograms directly. This is the same signal DatacardMaker uses to decide
        whether the categories are sliced, so the two cannot disagree.
        """
        return bool(self.get_config_data().get("binning"))

    def config_binning(self):
        """The datacard configuration's ``binning:`` block, filled in from
        BINNING_DEFAULTS for anything it does not declare."""
        declared = self.get_config_data().get("binning") or {}
        unknown = set(declared) - set(BINNING_PARAMS)
        if unknown:
            raise RuntimeError(
                f"{self.datacard_config_path()}: unknown key(s) in the 'binning' block: "
                f"{sorted(unknown)}. Known keys: {sorted(BINNING_PARAMS)}."
            )
        return {p: declared.get(p, BINNING_DEFAULTS[p]) for p in BINNING_PARAMS}

    def binning_params(self):
        """Effective binning: an explicitly given parameter wins, otherwise the
        configuration's value."""
        config = self.config_binning()
        return {
            p: (getattr(self, p) if getattr(self, p) is not None else config[p])
            for p in BINNING_PARAMS
        }

    def binning_diffs(self):
        """Binning parameters overridden away from what the configuration asks for."""
        config = self.config_binning()
        return {p: v for p, v in self.binning_params().items() if v != config[p]}

    def binning_parts(self):
        """Extra path level identifying a non-default binning, as a tuple to splat.

        law decides completeness from output existence alone, so without this a changed
        binning parameter silently reuses the previous shapes -- the task reports complete
        and the limits are computed from the old edges.

        Empty for the default binning, so nominal output paths (and the
        --multi-datacards globs that point at them) are exactly where they always were,
        and only scan variants gain a directory level. HistRebinTask writes
        binning_params.json into its output so a hash directory can be decoded later.
        """
        diffs = self.binning_diffs()
        if not diffs:
            return ()
        payload = ",".join(f"{k}={diffs[k]!r}" for k in sorted(diffs))
        return ("bin_" + hashlib.sha1(payload.encode()).hexdigest()[:10],)

    def store_parts(self):
        return (
            self.version,
            self.__class__.__name__,
            *self.binning_parts(),
            self.period,
        )

    def datacards_dir(self, era):
        """Local directory holding an era's datacards.

        The cards are produced on fs_default but must be real local files by the time
        combine sees them; ResonantLimitsTask mirrors them here and everything downstream
        (dhi's --multi-datacards globbing, the overlay plots) resolves against this path.
        """
        return os.path.join(
            self.ana_data_path(),
            self.version,
            "Datacards",
            *self.binning_parts(),
            era,
        )


class HistRebinTask(RebinningParamsTask, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 2.0)
    n_cpus = copy_param(HTCondorWorkflow.n_cpus, 1)

    # If set, discover bin edges from the combined statistics of this meta-era's
    # sub-eras (self.period must be one of them) -- for building the meta-era's
    # own combined limit. If unset, discover from self.period alone -- for a
    # standalone single-era limit. These need different edges (and thus different
    # output paths), since a single era's own statistics can't support the same
    # bins as the 4-era combination.
    meta_era = luigi.Parameter(default="")

    def get_discovery_eras(self):
        """Real eras summed together to discover DNN/HME bin edges: self.meta_era's
        sub-eras if set, otherwise just [self.period] (self.period's own stats)."""
        if self.meta_era:
            return self.get_era_groups()[self.meta_era]
        return [self.period]

    def discovery_hist_reqs(self):
        """{(era, variable): MergedHists} -- every input file this task reads."""
        return self.merged_hist_reqs(self.get_discovery_eras())

    def workflow_requires(self):
        return {
            f"MergedHists_{era}_{variable}": req
            for (era, variable), req in self.discovery_hist_reqs().items()
        }

    def requires(self):
        return list(self.discovery_hist_reqs().values())

    def create_branch_map(self):
        return {0: None}

    def output(self):
        # standalone: HistRebin/<period>; meta-era: HistRebin/<meta_era>/<period> --
        # these use different bin edges and can't share a path.
        #
        # Written to fs_default rather than the local data dir: the sliced shapes are
        # ~1 GB per version, which does not belong on an AFS work quota, and the only
        # consumer (CreateDatacardsTask) already stages its input down with localize().
        parts = [self.meta_era, self.period] if self.meta_era else [self.period]
        return self.output_dir_target(
            self.version, "HistRebin", *self.binning_parts(), *parts
        )

    def run(self):
        if not self.rebinning_enabled():
            # Without a `binning:` block the input is already 1D: process_category would
            # skip every histogram on GetDimension() != 2 and leave empty slice
            # directories behind, which law would then call a complete output.
            raise RuntimeError(
                f"{self.datacard_config_path()} declares no 'binning' block, so its input "
                "is not the 2D shapes this task rebins. CreateDatacardsTask reads such a "
                "configuration's merged histograms directly and does not require this task."
            )
        config = self.datacard_config_path()
        hist_rebin_py = os.path.join(
            self.ana_path(), "StatInference", "dc_make", "hist_rebin_2d.py"
        )
        targets = {key: req.output() for key, req in self.discovery_hist_reqs().items()}
        self.check_inputs(targets)

        # localize("w") on the output yields a local scratch dir and copies it to
        # fs_default on clean exit, so a crashed run leaves no half-written remote
        # output behind.
        binning = self.binning_params()
        with contextlib.ExitStack() as stack:
            base_dir_local = stack.enter_context(self.stage_inputs(targets))
            local_output = stack.enter_context(self.output().localize("w"))
            cmd = [
                # -u: stdout is a pipe here, so Python block-buffers it and the
                # per-mass progress lines are lost if the job is killed -- which is
                # exactly when they are needed to tell how far it got.
                "python3",
                "-u",
                hist_rebin_py,
                "--input",
                base_dir_local.abspath,
                "--output",
                local_output.abspath,
                "--config",
                config,
                "--era",
                self.period,
                "--discovery-eras",
                ",".join(self.get_discovery_eras()),
            ]
            # Built from the resolved values, not from self.<param>: those are None
            # wherever the datacard configuration is supplying the value.
            for name in BINNING_PARAMS:
                cmd += [f"--{name.replace('_', '-')}", str(binning[name])]
            ps_call(cmd, env=self.cmssw_env, verbose=1)

            # Decodes the bin_<hash> directory name back into the parameters that
            # produced it.
            with open(
                os.path.join(local_output.abspath, "binning_params.json"), "w"
            ) as f:
                json.dump(binning, f, indent=2, sort_keys=True)


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
                cmd += ["--n-dnn-slices", str(binning["n_dnn_slices"])]
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
            # The DNN slices span ~5 decades in yield (the low-significance slice holds
            # most of the background), so a linear axis hides everything but slice 0.
            "--log-y",
        ]
        if self.rebinning_enabled():
            # Same reason as create_datacards.py above: the panels are the sliced
            # categories, so they must follow the values HistRebinTask actually used.
            # Left unset otherwise, so the panels are the config's own categories.
            binning = self.binning_params()
            cmd += ["--n-dnn-slices", str(binning["n_dnn_slices"])]
            cmd += ["--category-pattern", binning["category_pattern"]]
        try:
            ps_call(cmd, verbose=1)
        except Exception as e:
            # The datacards are the product that matters; a plotting failure should be
            # loud but must not leave the task looking broken with valid cards on disk.
            print(
                f"WARNING: rebinned-shape plotting failed, datacards are unaffected: {e}"
            )


class ResonantLimitsTask(RebinningParamsTask):
    workflow = luigi.Parameter(default=law.parameter.NO_STR)

    def store_parts(self):
        return (
            self.version,
            self.__class__.__name__,
            *self.binning_parts(),
            "combined",
        )

    def get_eras(self):
        return self.get_top_level_eras()

    def _create_datacards_req(self, era, **kwargs):
        era_groups = self.get_era_groups()
        if era in era_groups:
            return CreateDatacardsTask.req(
                self, period=era_groups[era][0], meta_era=era, **kwargs
            )
        return CreateDatacardsTask.req(self, period=era, **kwargs)

    def requires(self):
        return [self._create_datacards_req(e, branches=()) for e in self.get_eras()]

    def output(self):
        return {
            "limits": self.local_target("limits.npz"),
            "datacards": law.LocalDirectoryTarget(self.datacards_dir("combined")),
        }

    def stage_datacards(self, era, remote_target):
        """Mirror an era's datacards from fs_default to a stable local path, returned.

        CreateDatacardsTask writes to fs_default, but the cards must be real local files
        by the time combine sees them: MergeResonantLimits shells out to combine, and
        PlotMultipleResonantLimits resolves --multi-datacards by globbing the filesystem
        outside law entirely. A localize() scratch dir will not do -- run() yields, which
        suspends and re-enters it, so the scratch dir is gone before dhi reads anything.

        The mirror path is the one CreateDatacardsTask used to write to directly, so
        existing --multi-datacards invocations keep working unchanged.
        """
        local_dir = self.datacards_dir(era)
        local_target = law.LocalDirectoryTarget(local_dir)
        # Re-staged on every run, including the re-entry after the yield below. The
        # copy is small (~30 MB) and unconditional refresh is what keeps the mirror
        # from going stale when CreateDatacardsTask reruns within the same version.
        if local_target.exists():
            local_target.remove()
        local_target.touch()
        remote_target.copy_to_local(local_target)
        return local_dir

    def run(self):
        datacards = []
        eras = self.get_eras()
        era_cards = {}

        for e in eras:
            create_dc_br0 = self._create_datacards_req(e, branch=0, branches=())
            output_dir = self.stage_datacards(e, create_dc_br0.output())
            cards = glob.glob(os.path.join(output_dir, "*.txt"))
            era_cards[e] = cards
            datacards.extend(cards)

        limits = yield MergeResonantLimits(
            version=self.version, datacards=tuple(datacards)
        )
        print(f"Merged limits: {limits}")

        self.output()["limits"].parent.touch()
        shutil.copy2(limits.path, self.output()["limits"].path)

        out_dc_dir = self.output()["datacards"]
        out_dc_dir.touch()

        masses = set()
        for e, cards in era_cards.items():
            for c in cards:
                m = re.search(r"_(\d+)\.txt$", c)
                if m:
                    masses.add(m.group(1))

        for mass in masses:
            combine_args = []
            for e in eras:
                for c in era_cards[e]:
                    if c.endswith(f"_{mass}.txt"):
                        combine_args.append(f"{e}={c}")
                        break

            if combine_args:
                cmd = ["combineCards.py"] + combine_args
                out_file = os.path.join(out_dc_dir.path, f"combined_{mass}.txt")
                with open(out_file, "w") as f:
                    subprocess.run(cmd, env=self.cmssw_env, stdout=f, check=True)


class PlotResonantLimitsTask(RebinningParamsTask):
    """The limit plots, declared in the datacard configuration.

    Each entry of the config's ``limit_plots`` block becomes one dhi
    PlotMultipleResonantLimits task: a list of datacard globs with the labels they should
    carry, optional external limit curves, and the plot styling. Globs are relative to the
    era's datacards directory and may use ``${ERA}``; the per-channel and per-category
    sub-directories they address come from DatacardMaker.writeDatacards.

    An entry with ``bands: true`` additionally gets one dhi PlotResonantLimits plot per
    curve -- the standard single-curve plot with the +-1/+-2 sigma bands, which the overlay
    cannot show. External limits and luminosity projections are separate curves rather than
    datacards, so they appear only on the overlay.

    The dhi tasks are invoked as subprocesses rather than yielded as dynamic
    dependencies. luigi round-trips a dynamic dependency through to_str_params() /
    from_str_params(), and dhi declares lumi_scale as FloatParameter(default=None), whose
    unset value serialises to "None" and then raises on float("None"). Going through the
    command line keeps dhi unmodified: the default never leaves the plot process.

    The finished plots are copied to fs_default alongside the rest of the chain's
    products. dhi still writes its own copy under inference/data/store, which is what makes
    an already-drawn plot cheap to re-request; the fs_default copy is the durable one.
    """

    workflow = luigi.Parameter(default=law.parameter.NO_STR)

    # Forwarded to the dhi task as "--remove-output 0,a,y": depth 0 drops the plot itself
    # and nothing below it, so no fit is recomputed.
    redraw = luigi.BoolParameter(
        default=False,
        significant=False,
        description="discard existing limit plots and draw them again",
    )

    def requires(self):
        # ResonantLimitsTask is what mirrors the cards to datacards_dir(), which is what
        # the globs below resolve against -- the plots cannot be built before it has run.
        return [ResonantLimitsTask.req(self)]

    def output(self):
        # fs_default, like HistRebinTask and CreateDatacardsTask. Holds one
        # <era>/ sub-directory of plots plus plots.json naming them.
        return self.output_dir_target(self.version, "LimitPlots", *self.binning_parts())

    def _resolve_config_path(self, path):
        return path if os.path.isabs(path) else os.path.join(self.ana_path(), path)

    def common_plot_params(self, entry, era):
        """The dhi parameters an entry's plots share, whatever kind they are.

        Validated against PlotMultipleResonantLimits, which is the subclass and so carries
        every parameter the single-curve PlotResonantLimits has plus its own -- one check
        covers both plot kinds, and the band builder drops what does not apply.
        """
        name = entry.get("name") or "limits"

        params = {
            # Distinguishes the plot files of entries that differ only in labelling;
            # not significant, so it does not move the output directory.
            "plot_postfix": name,
            # Per era, not per entry: the same entry is built for every top-level era and
            # each one carries a different luminosity. plot_params may still override it.
            "campaign": self.get_campaign(era),
        }
        # get_params() rather than get_param_names(), which drops the insignificant ones
        # -- that is most of the plot styling (x_min, campaign, ...).
        known = {
            param_name for param_name, _ in PlotMultipleResonantLimits.get_params()
        }
        for key, value in (entry.get("plot_params") or {}).items():
            if key not in known:
                raise RuntimeError(
                    f"limit_plots entry '{name}': '{key}' is not a "
                    "PlotMultipleResonantLimits parameter."
                )
            params[key] = value

        # dhi turns an unset campaign into None and simply draws no luminosity label, so
        # a missing entry would silently produce an unlabelled plot. Refuse instead.
        if not params.get("campaign"):
            raise RuntimeError(
                f"limit_plots entry '{name}': no campaign for era '{era}'. Add it to the "
                f"'campaigns' map in {self.datacard_config_path()} (a key of "
                "dhi.config.campaign_labels), or set it in this entry's plot_params."
            )
        if params["campaign"] not in campaign_labels:
            # dhi falls back to drawing the key itself, which is a usable escape hatch for
            # a one-off label but is also exactly what a typo looks like.
            print(
                f"WARNING: limit_plots entry '{name}': campaign '{params['campaign']}' is "
                "not in dhi.config.campaign_labels; it will be drawn verbatim as the "
                "luminosity label."
            )

        return params

    def build_plot_spec(self, entry, era, extra_external=()):
        """PlotMultipleResonantLimits parameters for a ``limit_plots`` entry and era.

        Kept as a plain dict so the same values can both construct the task (to learn
        where it will write) and be rendered into a command line (to make it write).

        ``extra_external`` holds external-limit files generated for this run (the
        luminosity projections). They come first, so a projection of our own curve is
        drawn before any fixed reference such as the Run 2 result.
        """
        sequences, names = self.datacard_sequences(entry, era)
        external = tuple(extra_external) + tuple(
            self._resolve_config_path(p) for p in entry.get("external_limits") or []
        )

        return dict(
            version=self.version,
            multi_datacards=tuple(sequences),
            datacard_names=tuple(names),
            external_limits=external,
            **self.common_plot_params(entry, era),
        )

    def build_band_specs(self, entry, era):
        """PlotResonantLimits parameters, one per curve, for an entry with ``bands: true``.

        The band plot is drawn from a single datacard set, so the multi-curve parameters
        (multi_datacards, datacard_names, external_limits, colors, markers) have no meaning
        here and are dropped by filtering against the parameters PlotResonantLimits
        actually declares. Everything else -- axes, xsec, campaign -- is shared with the
        overlay by construction, so the two plots of an entry cannot drift apart.
        """
        if not entry.get("bands"):
            return []

        sequences, names = self.datacard_sequences(entry, era)
        shared = self.common_plot_params(entry, era)
        known = {param_name for param_name, _ in PlotResonantLimits.get_params()}

        specs = []
        for sequence, label in zip(sequences, names):
            spec = {key: value for key, value in shared.items() if key in known}
            spec.update(
                version=self.version,
                datacards=tuple(sequence),
                # Entry name and curve label both. All three entries of the HH->bbWW
                # config draw "*.txt", so dhi would hash them into one output directory;
                # the postfix is what keeps their band plots from overwriting each other.
                # dhi's join_postfix strips whatever is not [a-zA-Z0-9._+-], so a label
                # like "e#mu" needs no sanitising here.
                plot_postfix=f"{shared['plot_postfix']}_{label}",
            )
            specs.append(spec)
        return specs

    def datacard_sequences(self, entry, era):
        """(datacard glob per curve, legend label per curve) for an entry, validated.

        Shared by the plot itself and by the luminosity projection, which has to scale the
        same cards the leading curve is drawn from.
        """
        name = entry.get("name") or "limits"
        base_dir = self.datacards_dir(era)

        sequences, names = [], []
        for card in entry["datacards"]:
            pattern = os.path.join(
                base_dir, Template(card["glob"]).safe_substitute(ERA=era)
            )
            if not glob.glob(pattern):
                raise RuntimeError(
                    f"limit_plots entry '{name}': no datacard matches '{pattern}'. "
                    f"Check the glob against the cards under {base_dir}."
                )
            label = card["name"]
            if "{" in label or "}" in label:
                # law brace-expands datacard_names, so "fb^{-1}" would arrive as "fb^-1".
                raise RuntimeError(
                    f"limit_plots entry '{name}': datacard name '{label}' contains braces, "
                    "which law brace-expands. Put the information in the campaign label instead."
                )
            sequences.append((pattern,))
            names.append(label)
        return sequences, names

    def limits_npz(self, sequence):
        """The MergeResonantLimits .npz for a datacard sequence, produced if absent.

        The path is taken from dhi's own task rather than assembled here: the store
        directory is a hash of the resolved datacard set, which is not something to
        reimplement. Passing the glob is equivalent to passing the resolved paths --
        dhi resolves it in modify_param_values before hashing.
        """
        target = MergeResonantLimits(
            version=self.version, datacards=tuple(sequence)
        ).output()
        if not os.path.exists(target.path):
            ps_call(
                [
                    "law",
                    "run",
                    "MergeResonantLimits",
                    "--version",
                    self.version,
                    "--datacards",
                    ",".join(sequence),
                ],
                cwd=self.ana_path(),
                verbose=1,
            )
        if not os.path.exists(target.path):
            raise RuntimeError(
                f"MergeResonantLimits produced no limits for {list(sequence)}; expected "
                f"{target.path}."
            )
        return target.path

    def make_lumi_projection(self, entry, era, projection, sequence, out_dir):
        """Write a dhi external-limits JSON scaling this entry's own limit curve to a
        different integrated luminosity, and return its path.

        Regenerated from the current limits on every run rather than read from a
        checked-in file: the projection is a function of the measured curve, so a stored
        copy silently goes stale the moment the binning, the datacards or the fits change.
        """
        name = entry.get("name") or "limits"
        label = projection["name"]
        campaign = self.get_campaign(era)
        # Defaults to the luminosity dhi prints on the plot for this campaign, so the
        # curve cannot be scaled from a different number than the one drawn beside it.
        lumi_now = projection.get("lumi_now", campaign_lumis.get(campaign))
        if not lumi_now:
            raise RuntimeError(
                f"limit_plots entry '{name}': lumi_projection '{label}' has no lumi_now "
                f"and campaign '{campaign}' is not in dhi.config.campaign_lumis."
            )

        out_path = os.path.join(
            out_dir, re.sub(r"[^\w.-]+", "_", f"{name}_{label}").strip("_") + ".json"
        )

        import numpy as np

        # The same purely statistical scaling dhi's own --lumi-scale applies
        # (plots/limits.py): limit(L_target) = limit(L_now) * sqrt(L_now / L_target).
        # It assumes the systematics scale away with the data, which they do not, so
        # this curve is optimistic -- an extrapolation of the statistical reach, not a
        # projected result. --lumi-scale itself is not usable here: it exists only on
        # plot_limit_scan, and it overwrites the measured curve rather than adding a
        # second one, whereas the point is to draw both side by side.
        #
        # Only "factor" is written, not pre-scaled limits: dhi's read_external_limits
        # (tasks/resonant.py) multiplies by factor * scale itself, so applying it here
        # too would double-count it.
        data = np.load(self.limits_npz(sequence), allow_pickle=True)["data"]
        entry_json = [
            {
                "name": label,
                "scan_parameter": "mhh",
                "factor": math.sqrt(lumi_now / projection["lumi_target"]),
                "limits": {repr(float(r["mhh"])): float(r["limit"]) for r in data},
            }
        ]
        with open(out_path, "w") as f:
            json.dump(entry_json, f, indent=2)
            f.write("\n")
        return out_path

    def plot_targets(self, task_cls, spec):
        """Where the dhi task will write, without scheduling it.

        Constructing the task is safe and is the only honest way to learn its output
        paths -- they are built from a hash of the datacard sequences plus several
        parameters, which is not something to reimplement here.
        """
        return law.util.flatten(task_cls(**spec).output())

    def plot_command(self, task_cls, spec):
        """`law run <dhi plot task> ...` for a spec."""
        cmd = ["law", "run", task_cls.__name__]
        for key, value in spec.items():
            flag = "--" + key.replace("_", "-")
            if isinstance(value, bool):
                # luigi bool parameters are set by presence, not by value
                if value:
                    cmd.append(flag)
            elif key == "multi_datacards":
                # colon between datacard sequences, comma within one
                cmd += [flag, ":".join(",".join(seq) for seq in value)]
            elif isinstance(value, (tuple, list)):
                cmd += [flag, ",".join(str(v) for v in value)]
            else:
                cmd += [flag, str(value)]
        if self.redraw:
            cmd += ["--remove-output", "0,a,y"]
        return cmd

    def draw_plot(self, task_cls, spec, name, era, dest_dir):
        """Run a dhi plot task, check it actually drew, copy the result into ``dest_dir``.

        Returns the basenames copied, for the plots.json manifest.
        """
        # cwd is pinned to the analysis root: dhi's resolve_datacards() takes a
        # different branch when the process happens to sit inside a configured
        # datacards_run2 directory.
        ps_call(self.plot_command(task_cls, spec), cwd=self.ana_path(), verbose=1)

        basenames = []
        for target in self.plot_targets(task_cls, spec):
            # luigi's retcode defaults return 0 even when a task fails, and no
            # [retcode] section overrides them here -- so a clean exit is not
            # evidence that anything was drawn. The file is.
            if not os.path.exists(target.path):
                raise RuntimeError(
                    f"limit_plots entry '{name}' ({era}): law exited cleanly "
                    f"but {target.path} was not written. See the output above "
                    "for the task that actually failed."
                )
            shutil.copy2(target.path, dest_dir)
            basenames.append(os.path.basename(target.path))
        return basenames

    def run(self):
        entries = self.get_config_data().get("limit_plots", [])
        if not entries:
            raise RuntimeError(
                f"{self.datacard_config_path()} declares no 'limit_plots' block, so there "
                "is nothing to plot."
            )

        with self.output().localize("w") as local_output:
            produced = []
            for era in self.get_eras_to_plot():
                dest_dir = os.path.join(local_output.abspath, era)
                os.makedirs(dest_dir, exist_ok=True)

                for entry in entries:
                    name = entry.get("name") or "limits"

                    # Projections are written into the output directory before the plot
                    # runs: they are --external-limits inputs, and they are worth keeping
                    # next to the plot as the record of what was actually drawn.
                    sequence = self.datacard_sequences(entry, era)[0][0]
                    projections = [
                        self.make_lumi_projection(entry, era, proj, sequence, dest_dir)
                        for proj in entry.get("lumi_projections") or []
                    ]

                    specs = [
                        (
                            PlotMultipleResonantLimits,
                            self.build_plot_spec(
                                entry, era, extra_external=projections
                            ),
                        )
                    ]
                    specs += [
                        (PlotResonantLimits, s)
                        for s in self.build_band_specs(entry, era)
                    ]

                    for task_cls, spec in specs:
                        produced.extend(
                            os.path.join(era, b)
                            for b in self.draw_plot(task_cls, spec, name, era, dest_dir)
                        )
                    produced.extend(
                        os.path.join(era, os.path.basename(p)) for p in projections
                    )

            with open(os.path.join(local_output.abspath, "plots.json"), "w") as f:
                json.dump(sorted(produced), f, indent=2)

    def get_eras_to_plot(self):
        return self.get_top_level_eras()


class ResonantLimitsAndHistPlotTask(RebinningParamsTask):
    workflow = luigi.Parameter(default=law.parameter.NO_STR)

    def get_eras(self):
        """Real, plottable periods only -- meta-era names are excluded since
        HistPlotTask requires an actual production period."""
        era_groups = self.get_era_groups()
        eras = self.get_config_data().get("eras", [self.period])
        return [e for e in eras if e not in era_groups]

    def requires(self):
        reqs = [ResonantLimitsTask.req(self)]
        for e in self.get_eras():
            reqs.append(HistPlotTask.req(self, period=e))
        return reqs

    def output(self):
        return self.local_target("dummy.txt")

    def run(self):
        self.output().touch()
