import argparse
import datetime
import json
import os
import re
import shutil
import subprocess
import sys
import yaml

file_dir = os.path.dirname(os.path.abspath(__file__))
pkg_dir = os.path.dirname(file_dir)
base_dir = os.path.dirname(pkg_dir)
pkg_dir_name = os.path.split(pkg_dir)[1]
if base_dir not in sys.path:
    sys.path.append(base_dir)
__package__ = pkg_dir_name

from FLAF.RunKit.run_tools import ps_call
from FLAF.RunKit.envToJson import get_cmsenv

cmssw_env = get_cmsenv(cmssw_path=os.getenv("FLAF_CMSSW_BASE")) 
for var in [ 'HOME', 'ANALYSIS_PATH', 'ANALYSIS_DATA_PATH', 'X509_USER_PROXY', 'CENTRAL_STORAGE',
            'ANALYSIS_BIG_DATA_PATH', 'FLAF_CMSSW_BASE', 'FLAF_CMSSW_ARCH' ]:
    if var in os.environ:
        cmssw_env[var] = os.environ[var]


def sh_call(cmd, error_message, verbose=0):
    if verbose > 0:
        print('>> {}'.format(cmd))
    returncode = subprocess.call([cmd], shell=True)
    if returncode != 0:
        raise RuntimeError(error_message)

def compare_binnings(b1, b2):
    if b2 is None: return True
    if b1['exp_limit'] != b2['exp_limit']: return b1['exp_limit'] < b2['exp_limit']
    if len(b1['bin_edges']) != len(b2['bin_edges']): return len(b1['bin_edges']) < len(b2['bin_edges'])
    for n in reversed(range(len(b1['bin_edges']))):
        if b1['bin_edges'][n] != b2['bin_edges'][n]: return b1['bin_edges'][n] < b2['bin_edges'][n]
    return False

def getBestBinning(log_file):
    with open(log_file, 'r') as f:
        binnings = json.loads('[' + ', '.join(f.readlines()) + ']')
    best_binning = None
    for binning in binnings:
        if compare_binnings(binning, best_binning):
            best_binning = binning
    return best_binning

def optimize_channel(channel, output, era, categories, max_n_bins, params, binning_suggestions, verbose):
    output_dir = os.path.join(output, channel+'_'+era)
    workers_dir = os.path.join(output_dir, 'workers')
    best_dir = os.path.join(output_dir, 'best')

    for d in [output, output_dir, workers_dir, best_dir]: 
        if not os.path.isdir(d):
            os.mkdir(d)

    suggested_binnings = {}

    for bs_file in binning_suggestions: 
        if os.path.isfile(bs_file):
            print(f"Loading suggested binnings from {bs_file}")
            with open(bs_file, 'r') as f:
                bs = json.load(f)
            if channel+'_'+era in bs: 
                for cat, cat_entry in bs[channel+'_'+era].items():
                    if cat not in suggested_binnings:
                        suggested_binnings[cat] = []
                    if type(cat_entry) == list:
                        if len(cat_entry) > 0:
                            if type(cat_entry[0]) == list:
                                for binning in cat_entry:
                                    suggested_binnings[cat].append(binning)
                            else:
                                suggested_binnings[cat].append(cat_entry)
                    elif type(cat_entry) == dict:
                        suggested_binnings[cat].append(cat_entry['bin_edges'])
                    else:
                        raise RuntimeError("Unknown format of suggested binning in '{}'.".format(bs_file))
        else:
            raise RuntimeError("Suggested_binnings file '{}' not found. Check cwd or fix bin_opt yaml file".format(bs_file))

    best_binnings_file = os.path.join(output_dir, 'best.json')
    if os.path.isfile(best_binnings_file):
        with open(best_binnings_file, 'r') as f:
            best_binnings = json.load(f)
    else:
        best_binnings = { channel+'_'+era: {} }

    first_cat_index = 0
    while first_cat_index < len(categories) and categories[first_cat_index][0] in best_binnings[channel+'_'+f"{era}"]:
        first_cat_index += 1

    for cat_index in range(first_cat_index, len(categories)):
        category, poi = categories[cat_index]
        print(f"Optimizing {channel} {category} for era {era}")
        for file in os.listdir(input_dir):
            if (channel in file or channel.lower() in file) and (category in file) and (era in file) and file.endswith('.txt'):
                input_card = f'{input_dir}/{file}'
            if (channel in file or channel.lower() in file) and (category in file) and (era in file) and file.endswith('.root'):
                input_shape = f'{input_dir}/{file}'
            else:
                continue

        cat_dir = os.path.join(output_dir, category)
        if not os.path.isdir(cat_dir):
            os.mkdir(cat_dir)

        cat_suggestions = os.path.join(cat_dir, 'to_try.json')
        if category in suggested_binnings:
            with open(cat_suggestions, 'w') as f:
                json.dump(suggested_binnings[category], f)

        cat_log = os.path.join(cat_dir, 'results.json')
        opt_cmd = f"python3 bin_opt/optimize_binning.py --input {input_card}  --shape-file {input_shape} --output {cat_dir} --workers-dir {workers_dir} --max_n_bins {max_n_bins} --poi {poi}"
        if params is not None: 
            opt_cmd += f" --params {params} " 
        for cat_idx in range(cat_index):
            cat = categories[cat_idx][0]
            for file in os.listdir(best_dir):
                if file.endswith('.txt') and cat in file:
                    other_cat_file = f"{best_dir}/{file}"

            if not os.path.isfile(other_cat_file):
                raise RuntimeError('Datacard "{}" for previous category not found.'.format(other_cat_file))
            opt_cmd += f' {other_cat_file} '

        ps_call(opt_cmd, shell=True, env=None, verbose=verbose)
        cat_best = getBestBinning(cat_log)
        if cat_best is None:
            raise RuntimeError("Unable to find best binning for {}".format(category))
        cat_best['poi'] = poi
        bin_edges = ', '.join([ str(edge) for edge in cat_best['bin_edges'] ])
        rebin_cmd = f'python3 bin_opt/rebinAndRunLimits.py --input {input_card} --output {best_dir} --bin-edges "{bin_edges}" --rebin-only '

        ps_call(rebin_cmd, shell=True, env=None, verbose=verbose)
        best_binnings[channel+'_'+era][category] = cat_best 
        with open(best_binnings_file, 'w') as f:
            f.write('{{\n\t"{}": {{\n'.format(channel+'_'+era))
            for cat_idx in range(0, cat_index + 1):
                cat = categories[cat_idx][0]
                f.write('\t\t "{}": '.format(cat))
                json.dump(best_binnings[channel+'_'+era][cat], f)
                if cat_idx < cat_index:
                    f.write(",")
                f.write("\n")
            f.write("\t}\n}\n")

    final_binning_file = output_dir + '.json'
    shutil.copy(best_binnings_file, final_binning_file)
    print("Binning for {} has been successfully optimised. The results can be found in {}" \
        .format(channel+'_'+era, final_binning_file)) 


if __name__ == "__main__":

    input_binning_opt_config = os.path.join(os.environ["ANALYSIS_PATH"], "StatInference", "bin_opt", "bin_optimization.yaml")
    with open(input_binning_opt_config, "r") as f:
        input_binning_opt_config_dict = yaml.safe_load(f)

    input_dir = input_binning_opt_config_dict["input"].get("directory", ".")
    channels = input_binning_opt_config_dict["input"].get("channels", [])
    era = str(input_binning_opt_config_dict["input"].get("era", ""))
    output = input_binning_opt_config_dict["output"].get("directory", "")
    max_n_bins = input_binning_opt_config_dict["input"].get("max_n_optimized_bins", 20)
    params = input_binning_opt_config_dict["input"].get("extra_parameters_optimize_binning", None)
    categories_poi = input_binning_opt_config_dict["input"].get("categories_ParameterOfInterest", "res2b:r")
    verbose = input_binning_opt_config_dict["output"].get("verbose", 0)
    binning_suggestions = input_binning_opt_config_dict["input"].get("binning_suggestions", [])

    categories = []
    for cat_entry in categories_poi.split(','):
        category, poi = cat_entry.split(':')
        categories.append([category, poi])
    
    for channel in channels:
        optimize_channel(channel, output, era, categories, max_n_bins, params, binning_suggestions, verbose)
        if verbose > 0:
            print(f"Finished optimizing channel: {channel}")