import os
import sys
import json
import argparse
import subprocess
import shutil
from FLAF.RunKit.run_tools import ps_call
from FLAF.RunKit.envToJson import get_cmsenv

cmssw_env = get_cmsenv(cmssw_path=os.getenv("FLAF_CMSSW_BASE")) #environment needs to be set up appropriately when GetLimits function is called to run law tasks such as UpperLimits or MergeResonantLimts
for var in [ 'HOME', 'ANALYSIS_PATH', 'ANALYSIS_DATA_PATH', 'X509_USER_PROXY', 'CENTRAL_STORAGE',
             'ANALYSIS_BIG_DATA_PATH', 'FLAF_CMSSW_BASE', 'FLAF_CMSSW_ARCH' ]:
    if var in os.environ:
        cmssw_env[var] = os.environ[var]

def compare_binnings(b1, b2):
    if b2 is None: return True
    if b1['exp_limit'] != b2['exp_limit']: return b1['exp_limit'] < b2['exp_limit']
    if len(b1['bin_edges']) != len(b2['bin_edges']): return len(b1['bin_edges']) < len(b2['bin_edges'])
    for n in reversed(range(len(b1['bin_edges']))):
        if b1['bin_edges'][n] != b2['bin_edges'][n]: return b1['bin_edges'][n] < b2['bin_edges'][n]
    return False

def getBestBinning(log_file):
    with open(log_file, 'r') as f:
        binnings = [json.loads(line) for line in f if line.strip()]
    best_binning = None
    for binning in binnings:
        if compare_binnings(binning, best_binning):
            best_binning = binning
    return best_binning


if __name__ == '__main__':
    file_dir = os.path.dirname(os.path.abspath(__file__))
    pkg_dir = os.path.dirname(file_dir)
    base_dir = os.path.dirname(pkg_dir)
    pkg_dir_name = os.path.split(pkg_dir)[1]
    if base_dir not in sys.path:
        sys.path.append(base_dir)
    __package__ = pkg_dir_name
    #These lines above ensure that the script can import sibling or parent modules/packages, even when executed as a standalone script, by adjusting the Python path and package context dynamically. This is especially useful in complex project structures or when running scripts from the command line.

    parser = argparse.ArgumentParser(description='Optimize binning for a given channel in an era')
    parser.add_argument('--input', required=True, type=str, help="input directory to datacards and shape files")
    parser.add_argument('--channel', required=True, type=str, help="channel_year, e.g. tauTau_2022")
    parser.add_argument('--output', required=True, type=str, help="output directory")
    parser.add_argument('--max-n-bins', required=False, type=int, default=20, help="maximum number of bins")
    parser.add_argument('--params', required=False, type=str, default=None, help="algorithm parameters in the format: param1=value1,param2=valu2,...")
    parser.add_argument('--verbose', required=False, type=int, default=1, help="verbosity level, 0: no output, 1: minimal output, 2: full output")
    parser.add_argument('--categories', required=False, type=str, default='res2b:r,res1b:r,boosted:r,classVBF:r_qqhh,classGGF:r,classttH:r_qqhh,classTT:r_qqhh,classDY:r_qqhh', help="comma separated list of categories and probability of improvement in the format --> category1:poi1,category2:poi2,...")
    parser.add_argument('--binning_suggestions', type=str, nargs='*', help="suggestions for binnings to try (e.g. best binnings from the previous round), json file")
    
    args = parser.parse_args()

    output_dir = os.path.join(args.output, args.channel)
    workers_dir = os.path.join(output_dir, 'workers')
    best_dir = os.path.join(output_dir, 'best')

    categories = [[category, poi] for category, poi in (cat_entry.split(':') for cat_entry in args.categories.split(','))]

    for dir in [args.output, output_dir, workers_dir, best_dir]:
        os.makedirs(dir, exist_ok=True)


    # def binning_suggestion_load_file(suggested_binnins={}, args.binning_suggestions, args.channel)
    suggested_binnings = {}

    for binning_suggestion_file in args.binning_suggestions:
        with open(binning_suggestion_file, 'r') as f:
            bs = json.load(f)
        if args.channel in bs:
            for cat, cat_entry in bs[args.channel].items():
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
                    raise RuntimeError(f"Unknown format of suggested binning in '{bs_file}' ")

    #def get_best_binnings(output_dir, args.channel)
    best_binnings_file = os.path.join(output_dir, 'best.json')
    if os.path.isfile(best_binnings_file):
        with open(best_binnings_file, 'r') as f:
            best_binnings = json.load(f)
    else:
        best_binnings = { args.channel: {}}


    #def (categories, best_binnings, args.channel)
    first_category_idx = 0
    while first_category_idx < len(categories) and categories[first_category_idx][0] in best_binnings[args.channel]:
        first_category_idx += 1

    for cat_index in range(first_category_idx, len(categories)):
        category, poi = categories[cat_index]
        print(f"Optimizing {args.channel} {category}")
        input_card = f'{args.input}/hh_{category}_{args.channel}_*.txt'
        cat_dir = os.path.join(output_dir, category)
        os.makedirs(cat_dir, exist_ok=True)

        category_suggestions = os.path.join(cat_dir, 'to_try.json')
        if category in suggested_binnings:
            with open(category_suggestions, 'w') as f:
                json.dump(suggested_binnings[category], f)

        category_log = os.path.join(cat_dir, 'results.json')

        # opt_cmd = f"python bin_opt/optimize_binning.py --input {input_card} --output {cat_dir} --workers-dir {workers_dir} --max-n-bins {args.max-n-bins} --poi {poi}"
        opt_cmd = [
            "python", "bin_opt/optimize_binning.py",
            "--input", input_card, "--output", cat_dir,
            "--workers-dir", workers_dir, "--max-n-bins", str(args.max_n_bins),
            "--poi", poi
        ]
        # opt_cmd += f" --params {args.params} " if args.params is not None else ""
        opt_cmd.append(f"--params={args.params}") if args.params is not None else None

        for cat_idx in range (cat_index):
            cat = categories[cat_idx][0]
            other_cat_file = f"{best_dir}/hh_{cat}_{args.channel}_*.txt"
            if not os.path.isfile(other_cat_file):
                raise RuntimeError(f"Datacard "{other_cat_file}" for category "{cat}" not found.")
            # opt_cmd += f' {other_cat_file} '
            opt_cmd.append(other_cat_file)
        
        # ps_call(opt_cmd, shell=True, env=cmssw_env, verbose=args.verbose)
        opt_proess = subprocess.Popen(opt_cmd, env=cmssw_env)

        cat_best = getBestBinning(category_log)
        if cat_best is None:
            raise RuntimeError(f"No best binning found for category {category} in {category_log}.")
        
        cat_best['poi'] = poi
        bin_edges = ', '.join([ str(edge) for edge in cat_best['bin_edges'] ])

        # ps_call(f"python bin_opt/rebinAndRunLimits.py --input {input_card} --output {best_dir} --bin-edges '{bin_edges}' --rebin-only ", shell=True, env=cmssw_env, verbose=args.verbose)
        rebin_cmd = [
            "python", "bin_opt/rebinAndRunLimits.py",
            "--input", input_card, "--output", best_dir,
            "--bin-edges", bin_edges, "--rebin-only"
        ]
        rebin_process = subprocess.Popen(rebin_cmd, env=cmssw_env)

        rebin_proc.wait()
        opt_proc.wait()

        best_binnings[args.channel][category] = cat_best
        with open(best_binnings_file, 'w') as f:
            f.write("{{\n\t'{}': {{\n".format(args.channel))
            for cat_idx in range(0, cat_index + 1):
                cat = categories[cat_idx][0]
                f.write(f"\t\t'{cat}': ")
                json.dump(best_binnings[args.channel][cat], f)
                if cat_idx < cat_index:
                    f.write(',')
                f.write("\n")
            f.write("\t}\n}\n")
    
    final_binning_file = output_dir + '.json'
    shutil.copy(best_binnings_file, final_binning_file)
    print(f"Binning for {args.channel} has been successfully optimized. The results can be found in {final_binning_file}")

