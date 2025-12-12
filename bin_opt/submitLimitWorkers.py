# import argparse
import os
# import subprocess
import yaml

from FLAF.RunKit.run_tools import ps_call
clean_env = {k: os.environ[k] for k in [
    'HOME', 'USER', 'LOGNAME', 'PATH', 'SHELL', 'ANALYSIS_SOFT_PATH', 'FLAF_CMSSW_BASE'
] if k in os.environ}


input_binning_opt_config = os.path.join(os.environ["ANALYSIS_PATH"], "StatInference", "bin_opt", "bin_optimization.yaml")
with open(input_binning_opt_config, "r") as f:
    input_binning_opt_config_dict = yaml.safe_load(f)

# parser = argparse.ArgumentParser(description='Submit rebinAndRunLimitsWorker jobs to HTCondor.')
# parser.add_argument('--output', required=True, type=str, help="output directory")
# parser.add_argument('--channel', required=True, type=str, help="channel")
# parser.add_argument('--suffix', required=False, type=str, default='', help="suffix to be added to the batch name")
# parser.add_argument('--n-workers', required=True, type=int, help="number of worker jobs to submit")
# parser.add_argument('--max-runtime', required=False, type=float, default=24, help="max runtime in hours")
# parser.add_argument('--verbose', required=False, type=int, default=1, help="verbosity")
# args = parser.parse_args()

output_dir = input_binning_opt_config_dict["output"]["directory"]
channels = input_binning_opt_config_dict["input"]["channels"]
era = input_binning_opt_config_dict["input"]["era"]
num_workers = input_binning_opt_config_dict["input"]["number_of_workers"]
max_runtime = input_binning_opt_config_dict["input"]["worker_condor_max_runtime_hours"]
verbose = input_binning_opt_config_dict["input"].get("verbose", 1)

for channel in channels:
    workers_output = os.path.join(output_dir, channel+'_'+str(era), "workers")
    print(f'workers_output {workers_output}')
    if not os.path.isdir(workers_output):
        raise RuntimeError("Workers output directory not found. Please, make sure that the server is running.")

    script_file = os.path.join(workers_output, 'run_worker.sh')

    with open(script_file, 'w') as f:
        f.write("python3 bin_opt/rebinAndRunLimitsWorker.py")
    os.chmod(script_file, 0o744)

    max_runtime_seconds = int(max_runtime * 60 * 60)

    condor_job = '''
    executable             = {0}
    log                    = {1}/worker.$(ClusterId).$(ProcId).log
    output                 = {1}/worker.$(ClusterId).$(ProcId).out
    error                  = {1}/worker.$(ClusterId).$(ProcId).err
    +MaxRuntime            = {2}
    queue {3}
    '''.format(script_file, workers_output, max_runtime_seconds, num_workers)

    sub_file = os.path.join(workers_output, 'worker.sub')
    with open(sub_file, 'w') as f:
        f.write(condor_job)
    
    ps_call("condor_submit -batch-name {} {}".format(channel, sub_file), env=clean_env, verbose=verbose, shell=True)

# workers_output = os.path.join(output_dir, args.channel, "workers")
# print(f'workers_output {workers_output}')
# if not os.path.isdir(workers_output):
#     raise RuntimeError("Workers output directory not found. Please, make sure that the server is running.")

# run_script = '''#!/bin/bash
# cd {}
# source env.sh hh
# python bin_opt/rebinAndRunLimitsWorker.py --output {} --verbose {}
# '''.format(os.getcwd(), workers_output, args.verbose)
# run_script = f'''
# python3 bin_opt/rebinAndRunLimitsWorker.py --output {workers_output} --verbose {args.verbose}
# '''
# script_file = os.path.join(workers_output, 'run_worker.sh')
# with open(script_file, 'w') as f:
#     f.write(run_script)
# os.chmod(script_file, 0o744)

# max_runtime_seconds = int(args.max_runtime * 60 * 60)

# condor_job = '''
# executable              = {0}
# log                     = {1}/worker.$(ClusterId).$(ProcId).log
# output                  = {1}/worker.$(ClusterId).$(ProcId).out
# error                   = {1}/worker.$(ClusterId).$(ProcId).err
# +MaxRuntime             = {2}
# queue {3}
# '''.format(script_file, workers_output, max_runtime_seconds, args.n_workers)

# sub_file = os.path.join(workers_output, 'worker.sub')
# with open(sub_file, 'w') as f:
#     f.write(condor_job)

# subprocess.call(['condor_submit -batch-name {} {}'.format(args.channel + args.suffix, sub_file)], shell=True)
