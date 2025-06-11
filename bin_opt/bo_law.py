import law
import luigi
import os
import json


from FLAF.RunKit.run_tools import ps_call #, natural_sort
from FLAF.RunKit.crabLaw import cond as kInit_cond, update_kinit_thread
from FLAF.run_tools.law_customizations import HTCondorWorkflow, copy_param, get_param_value #, Task
from FLAF.RunKit.envToJson import get_cmsenv

cmssw_env = get_cmsenv(cmssw_path=os.getenv("FLAF_CMSSW_BASE")) #environment needs to be set up appropriately when GetLimits function is called to run law tasks such as UpperLimits or MergeResonantLimts
for var in [ 'HOME', 'ANALYSIS_PATH', 'ANALYSIS_DATA_PATH', 'X509_USER_PROXY', 'CENTRAL_STORAGE',
             'ANALYSIS_BIG_DATA_PATH', 'FLAF_CMSSW_BASE', 'FLAF_CMSSW_ARCH' ]:
    if var in os.environ:
        cmssw_env[var] = os.environ[var]

class BinOptimizationTask(HTCondorWorkflow, law.LocalWorkflow): #(Task, HTCondorWorkflow, law.LocalWorkflow)
    max_runtime = copy_param(HTCondorWorkflow, 60.0)
    n_cpus = copy_param(HTCondorWorkflow, 4)

    def __init__(self, *args, **kwargs):
        super(BinOptimizationTask, self).__init__(*args, **kwargs)

    input_dir = luigi.Parameter()
    channel = luigi.Parameter()
    output_dir = luigi.Parameter()
    max_n_bins = luigi.IntParameter(default=20)
    params = luigi.Parameter(default="", significant=False)
    verbose = luigi.IntParameter(default=1)
    categories = luigi.Parameter(default="res2b:r,res1b:r,boosted:r,classVBF:r_qqhh,classGGF:r,classttH:r_qqhh,classTT:r_qqhh,classDY:r_qqhh")
    binning_suggestions = luigi.ListParameter(default=[])


    def create_branch_map(self):
        # Parse categories from string into a list of (category, poi) tuples
        cat_list = []
        for entry in self.categories.split(","):
            cat, poi = entry.split(":")
            cat_list.append({"category": cat, "poi": poi})
        return {i: cat_list[i] for i in range(len(cat_list))}

    def output(self):
        cat = self.branch_data["category"]
        return law.LocalFileTarget(os.path.join(self.output_dir, self.channel, cat, "best.json"))

    def run(self):
        cat = self.branch_data["category"]
        poi = self.branch_data["poi"]
        os.makedirs(os.path.dirname(self.output().path), exist_ok=True)
        cmd = [
            "python3", "bin_opt/optimize_distributor.py",
            "--input", self.input_dir,
            "--channel", self.channel,
            "--output", os.path.join(self.output_dir, self.channel, cat),
            "--max-n-bins", str(self.max_n_bins),
            "--categories", f"{cat}:{poi}",
            "--verbose", str(self.verbose)
        ]

        if self.params:
            if isinstance(self.params, dict):
                params_str = ",".join(f"{k}={v}" for k, v in self.params.items())
            elif isinstance(self.params, str):
                params_str = str(self.params)
            else:
                raise ValueError("params argument should be a dict `{param1=value1,param2=value2}` or a string `param1=value1,param2=value2.")
            cmd.extend( ["--params", params_str] )
        
        if self.binning_suggestions is not []:
            cmd.append("--binning_suggestions")
            cmd.extend(self.binning_suggestions)

        ps_call(" ".join(cmd), shell=True, env=cmssw_env, verbose=self.verbose)



# worker task to be submitted on one node
class RunRebinWorker(law.LocalWorkflow): #HTCondorWorkflow, 
    output_dir = luigi.Parameter()
    verbose = luigi.IntParameter(default=1)

    def output(self):
        # placeholder output to start the rebinning task
        return law.LocalFileTarget(os.path.join(self.output_dir, "worker.txt"))
    
    def run(self):
        cmd = [
            "python3", "bin_opt/rebinAndRunLimitsWorker.py",
            "--output", self.output_dir,
            "--verbose", str(self.verbose)
        ]
        ps_call(" ".join(cmd), shell=True)
        with self.output().open('w') as f:
            f.write("Worker started successfully.\n")


# optimizer task 
class RunOptimizeChannel(law.LocalWorkflow): # HTCondorWorkflow,
    input_dir = luigi.Parameter()
    output_dir = luigi.Parameter()
    channel = luigi.Parameter()
    verbose = luigi.IntParameter(default=1)

    def requires(self):
        # requires the worker task to have started
        return RunRebinWorker(output_dir=self.output_dir, verbose=self.verbose)

    def output(self):
        return law.LocalFileTarget(os.path.join(self.output_dir, f"{self.channel}_optimization.json"))
    
    def run(self):
        cmd = [
            "python3", "bin_opt/optimize_channel.py",
            "--input", self.input_dir,
            "--output", self.output_dir,
            "--channel", self.channel,
            "--verbose", str(self.verbose)
        ]
        ps_call(" ".join(cmd), shell=True)
        with self.output().open('w') as f:
            f.write("Channel optimization completed successfully.\n")

# Task to run both worker and optimizer tasks in parallel
class ParallelOptimizationWorkflow(law.WrapperTask):
    input_dir = luigi.Parameter()
    output_dir = luigi.Parameter()
    channel = luigi.Parameter()
    verbose = luigi.IntParameter(default=1)

    def requires(self):
        return [
            RunRebinWorker(output_dir=self.output_dir, verbose=self.verbose),
            RunOptimizeChannel(input_dir=self.input_dir, output_dir=self.output_dir, channel=self.channel, verbose=self.verbose)
        ]
#law run ParallelOptimizationWorkflow --input /path/to/input --output /path/to/output --channel my_channel --verbose 1