import numpy as np
import os
import uproot
import json
import sys
import matplotlib.pyplot as plt
from FLAF.RunKit.envToJson import get_cmsenv
from decimal import Decimal

# I had to do this because now the ps_call for MergeResonantLimits wouldn't work in cmsEnv environment
cmssw_env = get_cmsenv(cmssw_path=os.getenv("FLAF_CMSSW_BASE"))
for var in [
    "HOME",
    "ANALYSIS_PATH",
    "ANALYSIS_DATA_PATH",
    "X509_USER_PROXY",
    "CENTRAL_STORAGE",
    "ANALYSIS_BIG_DATA_PATH",
    "FLAF_CMSSW_BASE",
    "FLAF_CMSSW_ARCH",
]:
    if var in os.environ:
        cmssw_env[var] = os.environ[var]


def step_y(hists, signal_name, bkg_names, x_low, x_high, y_low, y_high):
    # Build lists of all new sums and x/y
    sum_signal = [0, 0, 0]
    sum_background = [0, 0, 0]
    new_y_low = [y_low, y_low, y_low - 1]
    new_y_high = [y_high, y_high + 1, y_high]

    sum_signal[0] = integrate_uproot(
        hists[signal_name], x_low, x_high, new_y_low[0], new_y_high[0]
    )[0]

    sum_signal[1] = integrate_uproot(
        hists[signal_name], x_low, x_high, new_y_low[1], new_y_high[1]
    )[0]
    sum_signal[2] = integrate_uproot(
        hists[signal_name], x_low, x_high, new_y_low[2], new_y_high[2]
    )[0]

    sum_background_current = 0

    for bkg_name in bkg_names:
        sum_background[0] += integrate_uproot(
            hists[bkg_name], x_low, x_high, new_y_low[0], new_y_high[0]
        )[0]

        sum_background[1] += integrate_uproot(
            hists[bkg_name], x_low, x_high, new_y_low[1], new_y_high[1]
        )[0]
        sum_background[2] += integrate_uproot(
            hists[bkg_name], x_low, x_high, new_y_low[2], new_y_high[2]
        )[0]

    significance = [0, 0, 0]
    significance[0] = sum_signal[0] / (sum_background[0] ** 0.5)
    significance[1] = sum_signal[1] / (sum_background[1] ** 0.5)
    significance[2] = sum_signal[2] / (sum_background[2] ** 0.5)

    print(f"Signficances s/sqrt(b) = {significance}")
    print(f"With total signals {sum_signal}")

    # Don't allow arg 0, but don't forget to add 1 to this new 'arg'
    best_sig_arg = np.argmax(significance[1:]) + 1

    print(f"New ylow/yhigh = {new_y_low[best_sig_arg]}/{new_y_high[best_sig_arg]}")

    return new_y_low[best_sig_arg], new_y_high[best_sig_arg]


def validate_bin(hists, bkg_names_main, bkg_names, x_low, x_high, y_low, y_high):
    for bkg_name in bkg_names_main:
        if integrate_uproot(hists[bkg_name], x_low, x_high, y_low, y_high)[0] <= 0:
            return False
    tot_background = 0
    for bkg_name in bkg_names:
        tot_background += integrate_uproot(
            hists[bkg_name], x_low, x_high, y_low, y_high
        )[0]
    if tot_background <= 0:
        return False
    print(f"Important backgrounds {bkg_names_main} had non-zero entries")
    print(f"Total background was {tot_background}")
    print(f"Passed validation for bin")

    return True


def integrate_uproot(hist, start_x, end_x, start_y=0, end_y=0):
    values = hist.values()
    errors = hist.errors()
    if len(values.shape) == 1:
        int_values = np.sum(values[start : end + 1])
        int_errors = np.sqrt(np.sum((errors[start : end + 1]) ** (2.0)))
    else:
        int_values = np.sum(values[start_x : end_x + 1, start_y : end_y + 1])
        int_errors = np.sqrt(
            np.sum((errors[start_x : end_x + 1, start_y : end_y + 1]) ** (2.0))
        )
    return [int_values, int_errors]


def convert_to_json(bins_dict, params):
    json_dict = []
    for key in bins_dict.keys():
        ch = key.split("_")[0]
        cat = "_".join(key.split("_")[1:])
        exists = False
        if not exists:
            entry = {
                "combined_bins": bins_dict[key],
                "channels": [ch],
                "categories": [cat],
            }
            for key, value in params.items():
                entry[key] = [value]
            json_dict.append(entry)

    return json_dict


# This will need to loop over files too
def optimize_binning(
    shapes_dir,
    filename,
    sig_string,
    params,
    mass_list,
    outdir,
    config,
    nBinsMaxX,
    nBinsMaxY,
    background_names,
    important_backgrounds,
    categories,
    channels,
    tag,
):
    os.makedirs(outdir, exist_ok=True)
    file = uproot.open(os.path.join(shapes_dir, filename))

    limit_history = {}
    limit_list = []

    best_limit = np.inf
    best_bins = None

    nBinsMaxX = nBinsMaxX
    nBinsMaxY = nBinsMaxY

    bkg_names = background_names

    priority_list = []
    bins_dict = {}
    for cat in categories:
        for ch in channels:
            priority_list.append((ch, cat))
            dict_key = f"{ch}_{cat}"
            # bins_dict[dict_key] = [0.0, 1.0]
            bins_dict[dict_key] = [
                {"y_bin": [np.float64(0.0), np.float64(2500.0)], "x_bins": [0, 1.0]}
            ]
            limit_history[dict_key] = {"nBins": [], "limits": []}

    for ch, cat in priority_list:
        ch_cat = f"{ch}/{cat}/"
        dict_key = f"{ch}_{cat}"
        hist_names = [
            key.split(";")[0][len(ch_cat) :]
            for key in file.keys()
            if key.startswith(ch_cat)
        ]
        hists = {hist_name: file[ch_cat + hist_name] for hist_name in hist_names}
        bkg_names_with_unc = []
        bkg_names_main = []  # Only the important backgrounds
        for bkg_name in bkg_names:
            for hist_name in hist_names:
                if hist_name.startswith(bkg_name):
                    bkg_names_with_unc.append(hist_name)
                    if bkg_name in important_backgrounds:
                        bkg_names_main.append(hist_name)

        print(f"Starting ch cat {ch_cat} with mass {mass_list}")

        signal_name = sig_string  # Maybe update to get from yaml later

        signal_values = hists[signal_name].values()
        xaxis, yaxis = hists[signal_name].axes
        total_signal = np.sum(signal_values)
        total_signal_y_slices = np.sum(signal_values, axis=0)
        x_max, y_max = signal_values.shape

        print(f"Total signal is {total_signal}")

        if total_signal == 0:
            print(f"This category doesn't have signal, skip")
            continue

        print(f"Total signal in y bins is {total_signal_y_slices}")

        signal_per_bin_y = total_signal / nBinsMaxY

        print(f"So we want {signal_per_bin_y} ({total_signal}/{nBinsMaxY})")

        y_upper_limit = np.max(np.nonzero(total_signal_y_slices))
        y_lower_limit = np.min(np.nonzero(total_signal_y_slices))

        print(f"Signal only has events in y range {y_lower_limit} to {y_upper_limit}")

        print(total_signal_y_slices[y_lower_limit : y_upper_limit + 1])

        y_edges = yaxis.edges()
        y_binning = set()
        y_binning_indexes = set()
        y_binning.add(y_edges[y_upper_limit])
        y_binning_indexes.add(int(y_upper_limit))
        current_upper_limit = y_upper_limit
        for low_y in range(y_upper_limit, y_lower_limit - 1, -1):
            if (
                np.sum(total_signal_y_slices[low_y : current_upper_limit + 1])
                >= signal_per_bin_y
            ):
                current_upper_limit = low_y
                y_binning.add(y_edges[low_y])
                y_binning_indexes.add(low_y)

        y_binning.add(y_edges[low_y])
        y_binning_indexes.add(low_y)
        y_binning = sorted(y_binning)
        y_binning_indexes = sorted(y_binning_indexes)

        print(y_binning)
        print(y_binning_indexes)

        final_binning = []

        x_edges = xaxis.edges()
        for y_index in range(len(y_binning) - 1):
            final_binning.append(
                {
                    "y_bin": [y_binning[y_index], y_binning[y_index + 1]],
                    "x_bins": [x_edges[0], x_edges[x_max]],
                }
            )

        # Have defined y binning, now do our 1d x binning inside these chunks

        for y_index in range(len(y_binning) - 1):
            print(
                f"Going to check the y range {y_binning[y_index]} to {y_binning[y_index+1]}"
            )
            print(
                f"Or in indexes, {y_binning_indexes[y_index]} to {y_binning_indexes[y_index+1]}"
            )

            # #Start the loop over nBins from 1 to nBinsMax
            signal_this_slice = np.sum(
                total_signal_y_slices[
                    y_binning_indexes[y_index] : y_binning_indexes[y_index + 1]
                ]
            )
            for nbins in range(6, nBinsMaxX):
                bin_content_goal = signal_this_slice / nbins

                print(
                    f"Channel {ch_cat} has a total of {signal_this_slice} entries in this Y, so with {nbins} bins our per-bin goal is {bin_content_goal}"
                )

                done = False
                nBins = len(hists[signal_name].values()) - 1
                right_edge = nBins
                x_binning = set()
                x_binning.add(x_edges[x_max])
                while not done:
                    for left_edge in range(right_edge, -1, -1):
                        current_total_signal = integrate_uproot(
                            hists[signal_name],
                            left_edge,
                            right_edge,
                            y_binning_indexes[y_index],
                            y_binning_indexes[y_index + 1] - 1,
                        )[0]
                        current_total_backgrounds = 0

                        for bkg_name in bkg_names:
                            current_total_backgrounds += integrate_uproot(
                                hists[bkg_name],
                                left_edge,
                                right_edge,
                                y_binning_indexes[y_index],
                                y_binning_indexes[y_index + 1] - 1,
                            )[0]

                        if current_total_signal >= bin_content_goal:
                            print(
                                f"Achieved bin content goal {current_total_signal} >= {bin_content_goal}, check if backgrounds are valid"
                            )
                            negative_flag = False
                            for bkg_name in bkg_names_main:
                                if (
                                    integrate_uproot(
                                        hists[bkg_name],
                                        left_edge,
                                        right_edge,
                                        y_binning_indexes[y_index],
                                        y_binning_indexes[y_index + 1] - 1,
                                    )[0]
                                    <= 0
                                ):
                                    # print(f"Failed, background {bkg_name} had zero or negative entries")
                                    negative_flag = True

                            if negative_flag:
                                if left_edge == 0:
                                    done = True
                                    # Don't add binning because we failed validation, just drop this low bin I guess
                                continue

                            # print("Passed the validation, append to binning and move on")
                            # I don't like having to round, but this is due to float64 != float32 (example 0.708)
                            x_binning.add(round(x_edges[left_edge], 3))
                            right_edge = left_edge - 1
                            done = True

                        if left_edge == 0:
                            done = True
                            print("Hit the left edge")
                            x_binning.add(x_edges[left_edge])
                            break

                x_binning = sorted(x_binning)
                final_binning[y_index]["x_bins"] = x_binning
                print(x_binning)
                print(y_binning[y_index], y_binning[y_index + 1])

                bins_dict[dict_key] = final_binning

                print("Final dict binnings are ", bins_dict)

                json_dict = convert_to_json(bins_dict, params)

                with open(
                    os.path.join(outdir, f"tmp_binning_{sig_string}.json"), "w"
                ) as f:
                    json.dump(json_dict, f, indent=4)

                # Create the datacards
                datacards_dir = os.path.join(outdir, f"tmp_datacards_{sig_string}")
                binning_file = os.path.join(outdir, f"tmp_binning_{sig_string}.json")

                ps_call(
                    f"python3 ../dc_make/create_datacards.py --input {shapes_dir} --output {datacards_dir} --config {config} --hist-bins ./{binning_file} --param_values {mass_list}",
                    shell=True,
                    env=cmssw_env,
                )
                # Find out where the limits will be saved
                output = ps_call(
                    f"law run MergeResonantLimits --version {tag} --datacards '{datacards_dir}/*.txt' --print-output 0,False",
                    shell=True,
                    catch_stdout=True,
                    split="\n",
                )[1]
                limit_file_name = output[-2]
                # Run the limit calculation
                ps_call(
                    f"law run MergeResonantLimits --version {tag} --datacards '{datacards_dir}/*.txt' --remove-output 3,a,y",
                    shell=True,
                )

                # Now load the output numpy array and get value [1] (mass, val, up1, down1, up2, down2)
                print(np.load(limit_file_name).files)
                data = np.load(limit_file_name)
                exp_limit = data[data.files[0]][0][1]
                print(f"Limit was {exp_limit}")

                limit_history[dict_key]["nBins"].append(nbins)
                limit_history[dict_key]["limits"].append(exp_limit)
                limit_list.append(exp_limit)

                # If adding a bin does not improve (limit stays the same) then just move on
                if (exp_limit <= best_limit) or (nbins == 1):
                    best_limit = exp_limit
                    best_bins = bins_dict.copy()

                else:
                    print(
                        f"Did not improve limits at nbins {nbins}, best limit was {best_limit}, and current was {exp_limit}"
                    )
                    # Need to replace the bin dict because it now has the 'worse' new binning
                    bins_dict = best_bins.copy()
                    break

    print("Finished the limit scan, lets look at the history")
    print(limit_history)

    print(limit_list)
    plt.plot(limit_list)
    plt.yscale("log")
    plot_filename = os.path.join(outdir, f"binopt_history_{sig_string}.pdf")
    plt.savefig(plot_filename)
    plt.close()


if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    pkg_dir = os.path.dirname(file_dir)
    base_dir = os.path.dirname(pkg_dir)
    pkg_dir_name = os.path.split(pkg_dir)[1]
    if base_dir not in sys.path:
        sys.path.append(base_dir)
    __package__ = pkg_dir_name

    from FLAF.RunKit.run_tools import ps_call
    import argparse

    parser = argparse.ArgumentParser(description="Optimize Binning.")
    parser.add_argument(
        "--input-shapes", required=True, type=str, help="input directory of shape files"
    )
    parser.add_argument(
        "--output",
        required=True,
        type=str,
        help="output directory for new binning, shapes, and datacards",
    )
    parser.add_argument("--config", required=True, type=str, help="configuration file")
    parser.add_argument(
        "--mass-list", required=True, type=str, help="list of mass values to scan over"
    )
    parser.add_argument(
        "--nBinsMaxY",
        required=False,
        type=int,
        default=5,
        help="Maximum number of bins to test",
    )
    parser.add_argument(
        "--nBinsMaxX",
        required=False,
        type=int,
        default=7,
        help="Maximum number of bins to test",
    )
    parser.add_argument(
        "--sig-name-pattern",
        required=False,
        type=str,
        default="ggRadion_HH_bbWW_M{}",
        help="Name format of signal in the shape files",
    )
    parser.add_argument(
        "--file-name-pattern",
        required=False,
        type=str,
        default="Run3_2022/shape_m{m}.root",
        help="Name format of input shape files",
    )
    parser.add_argument(
        "--background-names",
        required=False,
        type=str,
        default="TT",
        help="List of all background processes in the shapes file",
    )
    parser.add_argument(
        "--important-backgrounds",
        required=False,
        type=str,
        default="TT",
        help="List of background processes to consider for individual statistics",
    )
    parser.add_argument(
        "--categories",
        required=False,
        type=str,
        default="res2b,boost,res1b",
        help="List of jet categories",
    )
    parser.add_argument(
        "--channels",
        required=False,
        type=str,
        default="eE,eMu,muMu",
        help="List of lepton channels",
    )
    parser.add_argument(
        "--tag", required=False, type=str, default="dev", help="Name of combine job"
    )

    args = parser.parse_args()

    mass_list = [int(x) for x in args.mass_list.split(",")]
    shapes_dir = args.input_shapes
    outdir = args.output
    config = args.config
    nBinsMaxY = args.nBinsMaxY
    nBinsMaxX = args.nBinsMaxX
    signal_name_pattern = args.sig_name_pattern
    file_name_pattern = args.file_name_pattern
    background_names = [x for x in args.background_names.split(",")]
    important_backgrounds = [x for x in args.important_backgrounds.split(",")]
    categories = [x for x in args.categories.split(",")]
    channels = [x for x in args.channels.split(",")]
    tag = args.tag

    for mass in mass_list:
        hists = optimize_binning(
            shapes_dir,
            file_name_pattern.format(m=mass),
            signal_name_pattern.format(m=mass),
            {"MX": mass},
            mass,
            outdir,
            config,
            nBinsMaxX,
            nBinsMaxY,
            background_names,
            important_backgrounds,
            categories,
            channels,
            tag,
        )
        # break
