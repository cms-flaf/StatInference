import json
import os
import shutil
import time
import uuid

from rebinAndRunLimits import GetLimits

import argparse

parser = argparse.ArgumentParser(description='Rebin histogram and run expected limits.')
parser.add_argument('--output', required=True, type=str, help="output directory")
parser.add_argument('--verbose', required=False, type=int, default=1, help="verbosity")
args = parser.parse_args()

worker_dir = os.path.join(args.output, uuid.uuid4().hex)
os.mkdir(worker_dir)
task_file = os.path.join(worker_dir, 'task.txt')
result_file = os.path.join(worker_dir, 'result.txt')
result_file_tmp = os.path.join(worker_dir, '.result.txt')

if args.verbose > 0:
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
        if args.verbose > 0:
            print('Bin edges: [ {} ]'.format(', '.join([ str(b) for b in bin_edges ])))
        limit = GetLimits(params['input_datacard'], worker_dir, bin_edges, params['poi'], verbose=1,
                          other_datacards=params['other_datacards'])
        print('back in rebinAndRunLimitsWorker after running limits')
        if args.verbose > 0:
            print('Expected 95% CL limit = {}'.format(limit))
        result = {
            'input_datacard': params['input_datacard'],
            'bin_edges': bin_edges,
            'exp_limit': limit
        }
        with open(result_file_tmp, 'w') as f:
            json.dump(result, f)
        shutil.move(result_file_tmp, result_file)
