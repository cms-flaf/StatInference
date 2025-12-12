import json
import os
import shutil
import time
import uuid

import yaml
input_binning_opt_config = os.path.join(os.environ["ANALYSIS_PATH"], "StatInference", "bin_opt", "bin_optimization.yaml")
with open(input_binning_opt_config, "r") as f:
    input_binning_opt_config_dict = yaml.safe_load(f)

from rebinAndRunLimits import GetLimits

for channel in input_binning_opt_config_dict["input"]["channels"]:
    if os.path.isfile(os.path.join(input_binning_opt_config_dict["output"]["directory"],
                                   f'{channel}_{str(input_binning_opt_config_dict["input"].get("era", ""))}.json')):
                                   print( f"Results for channel {channel} already exist.")
    else:
        worker_dir = os.path.join(input_binning_opt_config_dict["output"]["directory"],
        f'{channel}_{str(input_binning_opt_config_dict["input"].get("era", ""))}',
                                            "workers",
                                            uuid.uuid4().hex)

        os.makedirs(worker_dir, exist_ok=True)
        task_file = os.path.join(worker_dir, 'task.txt')
        result_file = os.path.join(worker_dir, 'result.txt')
        result_file_tmp = os.path.join(worker_dir, '.result.txt')

        verbose = input_binning_opt_config_dict.get("verbose", 1)

        if verbose > 0:
            print('Worker dir: {}'.format(worker_dir))

        while True:
            time.sleep(1)
            if os.path.isfile(task_file) and not os.path.isfile(result_file):
                try:
                    with open(task_file, 'r') as f:
                        params = json.load(f)
                except IOError:
                    continue
                bin_edges = params['bin_edges']
                if verbose > 0:
                    print('Bin edges: [ {} ]'.format(', '.join([ str(b) for b in bin_edges ])))
                limit = GetLimits(params['input_datacard'], worker_dir, bin_edges, params['poi'], verbose=1,
                                    other_datacards=params['other_datacards'])

                if verbose > 0:
                    print('Expected 95% CL limit = {}'.format(limit))
                result = {
                    'input_datacard': params['input_datacard'],
                    'bin_edges': bin_edges,
                    'exp_limit': limit
                }
                with open(result_file_tmp, 'w') as f:
                    json.dump(result, f)
                shutil.move(result_file_tmp, result_file)
