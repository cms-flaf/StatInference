import argparse
import datetime
import json
import os
import re
import shutil
import subprocess
import sys
import yaml
import math

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

def _total_bins_in_multicategory_entry(entry):
    categories = entry.get("categories", [])
    return sum(len(cat.get("bin_edges", [])) - 1 for cat in categories)

def getBestMultiCategoryResult(log_file):
    best = None
    with open(log_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            entry = json.loads(line)
            if best is None:
                best = entry
                continue
            if entry['exp_limit'] != best['exp_limit']:
                if entry['exp_limit'] < best['exp_limit']:
                    best = entry
            else:
                if _total_bins_in_multicategory_entry(entry) < _total_bins_in_multicategory_entry(best):
                    best = entry
    return best

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

def inputs_for_category(input_dir, channel, category, era):
    cards = []
    shapes = []
    for file in os.listdir(input_dir):
        if not ((channel in file or channel.lower() in file) and (category in file) and (era in file)):
            continue
        if file.endswith(".txt"):
            cards.append(os.path.join(input_dir, file))
        elif file.endswith(".root"):
            shapes.append(os.path.join(input_dir, file))
    if len(cards) != 1:
        raise RuntimeError(
            f"Expected exactly 1 datacard for channel={channel}, category={category}, era={era}, found {len(cards)}: {cards}"
        )
    if len(shapes) != 1:
        raise RuntimeError(
            f"Expected exactly 1 shape file for channel={channel}, category={category}, era={era}, found {len(shapes)}: {shapes}"
        )
    return os.path.abspath(cards[0]), os.path.abspath(shapes[0])

def optimize_channel(channel, output, era, categories, max_n_bins, params, binning_suggestions, verbose):
    output_dir = os.path.join(output, channel+'_'+era)
    workers_dir = os.path.join(output_dir, 'workers')
    best_dir = os.path.join(output_dir, 'best')

    for d in [output, output_dir, workers_dir, best_dir]: 
        if not os.path.isdir(d):
            os.makedirs(d)

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

    if input_binning_opt_config_dict["input"]["multi_category_optimization"]:
        remaining = categories[first_cat_index:]
        if len(remaining) == 0:
            print(f"Nothing to do: all categories already optimized for {channel}_{era}")
            return

        pois = sorted(set([poi for (_cat, poi) in remaining]))
        if len(pois) != 1:
            raise RuntimeError(f"Simultaneous multi-category optimization requires a single shared POI across the bundle. Found POIs: {pois}")
        poi = pois[0]

        bundle_categories = [cat for (cat, _poi) in remaining]
        datacards = []
        shape_files = []
        cat_inputs = {}
        for category in bundle_categories:
            card, shape = inputs_for_category(input_dir, channel, category, era)
            datacards.append(card)
            shape_files.append(shape)
            cat_inputs[category] = (card, shape)

        bundle_name = "multicat_" + "_".join(bundle_categories)
        bundle_dir = os.path.join(output_dir, bundle_name)
        os.makedirs(bundle_dir, exist_ok=True)

        input_arg = ",".join(datacards)
        shape_arg = ",".join(shape_files)
        opt_cmd = (
            f"python3 bin_opt/multi_optimize_binning.py "
            f"--input {input_arg} "
            f"--shape-file {shape_arg} "
            f"--output {bundle_dir} "
            f"--workers-dir {workers_dir} "
            f"--max_n_bins {max_n_bins} "
            f"--poi {poi}"
        )
        if params is not None:
            opt_cmd += f" --params {params} "

        ps_call(opt_cmd, shell=True, env=None, verbose=verbose)

        multicat_log = os.path.join(bundle_dir, "results_multibinning.json")
        if not os.path.isfile(multicat_log):
            raise RuntimeError(f"Multi-category log file not found: {multicat_log}")

        best_entry = getBestMultiCategoryResult(multicat_log)
        if best_entry is None:
            raise RuntimeError(f"Unable to find best multi-category result in {multicat_log}")

        combined_exp_limit = best_entry["exp_limit"]

        cats_sorted = sorted(best_entry["categories"], key=lambda x: int(x["cat_id"]))
        if len(cats_sorted) != len(bundle_categories):
            raise RuntimeError(
                "Mismatch between optimized categories and result categories: "
                f"{len(bundle_categories)} vs {len(cats_sorted)}"
            )

        for idx, (category, _poi) in enumerate(remaining):
            cat_result = cats_sorted[idx]
            bin_edges = cat_result["bin_edges"]

            input_card, _input_shape = cat_inputs[category]

            cat_best = {
                "bin_edges": bin_edges,
                "exp_limit": combined_exp_limit,  # combined objective
                "poi": poi,
                "mode": "simultaneous_multicat",
                "multicat_categories": bundle_categories,
            }
            best_binnings[channel + "_" + era][category] = cat_best

            bin_edges_str = ", ".join([str(edge) for edge in bin_edges])
            rebin_cmd = (
                f'python3 bin_opt/rebinAndRunLimits.py '
                f'--input {input_card} '
                f'--output {best_dir} '
                f'--bin-edges "{bin_edges_str}" '
                f'--rebin-only '
            )
            ps_call(rebin_cmd, shell=True, env=None, verbose=verbose)

        best_binnings_file = os.path.join(output_dir, "best.json")
        with open(best_binnings_file, "w") as f:
            f.write('{{\n\t"{}": {{\n'.format(channel + "_" + era))

            first_written = True
            for i, (category, _poi) in enumerate(categories):
                if category not in best_binnings[channel + "_" + era]:
                    continue
                if not first_written:
                    f.write(",\n")
                first_written = False
                f.write('\t\t "{}": '.format(category))
                json.dump(best_binnings[channel + "_" + era][category], f)
            if first_written:
                f.write("\t}\n}\n")
            else:
                f.write("\n\t}\n}\n")

        final_binning_file = output_dir + ".json"
        shutil.copy(best_binnings_file, final_binning_file)
        print(
            "Binning for {} has been successfully optimised (simultaneous multi-category). "
            "The results can be found in {}".format(channel + "_" + era, final_binning_file)
        )
        return


    for cat_index in range(first_cat_index, len(categories)):
        category, poi = categories[cat_index]
        print(f"Optimizing {channel} {category} for era {era}")
        datacards = []
        shape_files = []
        for file in os.listdir(input_dir):
            if (channel in file or channel.lower() in file) and (category in file) and (era in file) and file.endswith('.txt'):
                input_card = f'{input_dir}/{file}'
                datacards.append(f'{input_dir}/{file}')
            if (channel in file or channel.lower() in file) and (category in file) and (era in file) and file.endswith('.root'):
                input_shape = f'{input_dir}/{file}'
                shape_files.append(f'{input_dir}/{file}')
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
        if input_binning_opt_config_dict["input"]["multi_category_optimization"]:
            input_arg = ",".join([os.path.abspath(p) for p in datacards])
            shape_arg = ",".join([os.path.abspath(p) for p in shape_files])
            if input_arg == "" or shape_arg == "":
                raise RuntimeError("No datacards or shape files found in input directory")
            opt_cmd = f"python3 bin_opt/multi_optimize_binning.py --input {input_arg} --shape-file {shape_arg} --output {cat_dir} --workers-dir {workers_dir} --max_n_bins {max_n_bins} --poi {poi}"

        else:
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
    # channels = input_binning_opt_config_dict["input"].get("channels", [])
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
    
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--channel', default="tauTau", type=str, help="channel to optimize")
    args = parser.parse_args()
    channel = args.channel
    
    optimize_channel(channel, output, era, categories, max_n_bins, params, binning_suggestions, verbose)
    if verbose > 0:
        print(f"Finished optimizing channel: {channel}")
