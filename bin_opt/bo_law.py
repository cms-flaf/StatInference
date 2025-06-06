import law
import os


from FLAF.RunKit.run_tools import ps_call #, natural_sort
from FLAF.RunKit.crabLaw import cond as kInit_cond, update_kinit_thread
from FLAF.run_tools.law_customizations import Task, HTCondorWorkflow, copy_param, get_param_value
from FLAF.Common.Utilities import SerializeObjectToString
from FLAF.RunKit.envToJson import get_cmsenv

cmssw_env = get_cmsenv(cmssw_path=os.getenv("FLAF_CMSSW_BASE")) #environment needs to be set up appropriately when GetLimits function is called to run law tasks such as UpperLimits or MergeResonantLimts
for var in [ 'HOME', 'ANALYSIS_PATH', 'ANALYSIS_DATA_PATH', 'X509_USER_PROXY', 'CENTRAL_STORAGE',
             'ANALYSIS_BIG_DATA_PATH', 'FLAF_CMSSW_BASE', 'FLAF_CMSSW_ARCH' ]:
    if var in os.environ:
        cmssw_env[var] = os.environ[var]

class BinOptimizationTask(Task, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow, 60.0)
    n_cpus = copy_param(HTCondorWorkflow, 4)

    input_dir = luigi.Parameter()
    channel = luigi.Parameter()
    output_dir = luigi.Parameter()
    max_n_bins = luigi.IntParameter(default=20)
    params = luigi.Parameter(default=None, significant=False)
    verbose = luigi.IntParameter(default=1)
    categories = luigi.ListParameter(default="res2b:r,res1b:r,boosted:r,classVBF:r_qqhh,classGGF:r,classttH:r_qqhh,classTT:r_qqhh,classDY:r_qqhh")
    binning_suggestions = luigi.ListParameter(default=[])

    #def workflow_requires(self):
    #    return
    
    #def requires(self):
    #    return

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
            "python", "bin_opt/optimize_distributor.py",
            "--input", self.input_dir,
            "--channel", self.channel,
            "--output", os.path.join(self.output_dir, self.channel, cat),
            "--max-n-bins", str(self.max_n_bins),
            "--params", SerializeObjectToString(self.params), #for extra algorithm parameters in the future
            "--categories", f"{cat}:{poi}",
            "--verbose", str(self.verbose)
        ]

        ps_call(" ".join(cmd), shell=True, env=cmssw_env, verbose=self.verbose)

