import contextlib
import json
import law
import luigi
import os

from FLAF.RunKit.run_tools import ps_call
from FLAF.run_tools.law_customizations import HTCondorWorkflow, copy_param

from .RebinningParamsTask import BINNING_PARAMS, RebinningParamsTask


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
        """Real eras summed together to discover the bin edges: self.meta_era's
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
