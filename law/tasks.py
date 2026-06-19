import law
import luigi
import os

from FLAF.RunKit.run_tools import ps_call
from FLAF.run_tools.law_customizations import (
    Task,
    HTCondorWorkflow,
    copy_param,
)
from FLAF.Analysis.tasks import HistMergerTask, HistPlotTask
from dhi.tasks.resonant import MergeResonantLimits


class CreateDatacardsTask(Task, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 2.0)
    n_cpus = copy_param(HTCondorWorkflow.n_cpus, 1)

    def workflow_requires(self):
        return { "HistMerger": HistMergerTask.req(self, branches=()) }

    def requires(self):
        merge_map = HistMergerTask.req(self, branch=-1, branches=()).create_branch_map()

        return [
            HistMergerTask.req(self, branch=br, branches=(br,))
            for br in merge_map.keys()
        ]

    def create_branch_map(self):
        return { 0: None }

    def output(self):
        path = os.path.join(self.ana_data_path(), self.version, "Datacards", self.period)
        return law.LocalDirectoryTarget(path)

    def run(self):
        statInf_entry = self.global_params["StatInference"]
        config = os.path.join(self.ana_path(), statInf_entry["config"])
        hist_bins = os.path.join(self.ana_path(), statInf_entry["hist_bins"])
        param_values = statInf_entry.get("param_values", [])
        create_datacards_py = os.path.join(self.ana_path(), "StatInference", "dc_make", "create_datacards.py")
        base_input_dir_remote = self.input()[0].parent.parent.parent
        with base_input_dir_remote.localize("r") as base_dir_local:
            cmd = [
                "python3",
                create_datacards_py,
                "--input",
                base_dir_local.abspath,
                "--output",
                self.output().abspath,
                "--config",
                config,
                "--hist-bins",
                hist_bins,
                "--eras",
                self.period,
            ]
            if len(param_values) > 0:
                param_values_str = ",".join(str(v) for v in param_values)
                cmd += ["--param_values", param_values_str]
            ps_call(cmd, env=self.cmssw_env, verbose=1)


class ResonantLimitsTask(Task):
    workflow = luigi.Parameter(default=law.parameter.NO_STR)

    def get_eras(self):
        statInf_entry = self.global_params["StatInference"]
        config = os.path.join(self.ana_path(), statInf_entry["config"])
        import yaml
        with open(config, "r") as f:
            data = yaml.safe_load(f)
        return data.get("eras", [self.period])

    def requires(self):
        return [ CreateDatacardsTask.req(self, period=e, branches=()) for e in self.get_eras() ]

    def output(self):
        return self.local_target("dummy.txt")

    def run(self):
        datacards = []
        for e in self.get_eras():
            create_dc_br0 = CreateDatacardsTask.req(self, period=e, branch=0, branches=())
            output_dir = create_dc_br0.output().abspath
            datacards.append(os.path.join(output_dir, "*.txt"))

        limits = yield MergeResonantLimits(version=self.version, datacards=tuple(datacards))
        print(f"Merged limits: {limits}")
        self.output().touch()


class ResonantLimitsAndHistPlotTask(Task):
    workflow = luigi.Parameter(default=law.parameter.NO_STR)
    
    def get_eras(self):
        statInf_entry = self.global_params["StatInference"]
        config = os.path.join(self.ana_path(), statInf_entry["config"])
        import yaml
        with open(config, "r") as f:
            data = yaml.safe_load(f)
        return data.get("eras", [self.period])

    def requires(self):
        reqs = [ ResonantLimitsTask.req(self) ]
        for e in self.get_eras():
            reqs.append(HistPlotTask.req(self, period=e))
        return reqs

    def output(self):
        return self.local_target("dummy.txt")

    def run(self):
        self.output().touch()
