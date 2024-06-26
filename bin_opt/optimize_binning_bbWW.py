import numpy as np
import os
import uproot
import json
import sys
import matplotlib.pyplot as plt

def integrate_uproot(hist, start, end):
    int_values = np.sum(hist.values()[start:end+1])
    int_errors = np.sqrt(np.sum((hist.errors()[start:end+1])**(2.0)))
    return [int_values, int_errors]

def convert_to_json(bins_dict, params):
    json_dict = []
    for key in bins_dict.keys():
        ch = key.split("_")[0]
        cat = key.split("_")[1]
        exists = False
        #for entry in json_dict:
        #    if entry['bins'] == bins_dict[key]:
        #        entry['channels'].append(ch)
        #        entry['categories'].append(cat)
        #        exists = True
        if not exists:
            entry = {
                    'bins': bins_dict[key],
                    'channels': [ch],
                    'categories': [cat],
                }
            for key, value in params.items():
                entry[key] = [ value ]
            json_dict.append(entry)
                
            
    print(json_dict)
    return json_dict



#This will need to loop over files too
def optimize_binning(filename, params, sig_string, mass_list):
    file = uproot.open(filename)

    limit_history = {}
    limit_list = []

    best_limit = np.inf
    best_bins = None

    nBinsMax = 20

    bkg_names = ['TT', 'DY', 'ST', 'VV']

    priority_list = []
    bins_dict = {}
    for cat in [ "res2b", "boost", "res1b" ]:
        for ch in [ "mumu", "emu", "ee" ]:
            priority_list.append((ch, cat))
            dict_key = f'{ch}_{cat}'
            bins_dict[dict_key] = [0.0, 1.0]
            limit_history[dict_key] = {'nBins': [], 'limits': []}

    for ch, cat in priority_list:
        ch_cat = f'{ch}/{cat}/'
        dict_key = f'{ch}_{cat}'
        hist_names = [ key.split(';')[0][len(ch_cat):] for key in file.keys() if key.startswith(ch_cat) ]
        hists = { hist_name : file[ch_cat+hist_name] for hist_name in hist_names }

        print(f"Starting ch cat {ch_cat}")
        
        #signal_name = 'ggRadion_HH_bbWW_M300' #Update to get from yaml later
        signal_name = sig_string
        #important_backgrounds = ['TT', 'DY', 'ST']
        important_backgrounds = ['TT']


        tmp_signal = np.sum(hists[signal_name].values())

        #Start the loop over nBins from 1 to nBinsMax
        for nbins in range(1, nBinsMax):
            bin_content_goal = tmp_signal/nbins

            print(f"Channel {ch_cat} has a total of {tmp_signal} entries, so with {nbins} bins our per-bin goal is {bin_content_goal}")

            custom_binning = [1.0]
            custom_totals = []
            done = False
            nBins = len(hists[signal_name].values())-1
            right_edge = nBins
            big_bin_counter = 1
            while not done:
                for left_edge in range(right_edge, -1, -1):
                    tot_integral_signal = integrate_uproot(hists[signal_name], left_edge, nBins)[0]
                    curr_integral_signal = integrate_uproot(hists[signal_name], left_edge, right_edge)[0]
                    integral_bkgs = np.zeros(len(important_backgrounds))
                    integral_bkgs_unc = np.zeros(len(important_backgrounds))
                    tot_bkgs = 0
                    tot_bkgs_unc = 0
                    for i, bkg_name in enumerate(important_backgrounds):
                        int_bkg = integrate_uproot(hists[bkg_name], left_edge, right_edge)
                        if bkg_name == 'TT': #Only use TT for the later error check for now
                            integral_bkgs[i] = int_bkg[0]
                            integral_bkgs_unc[i] = (int_bkg[1])

                        tot_bkgs += int_bkg[i]
                        tot_bkgs_unc += (integral_bkgs_unc[i]**(2.0))

                    tot_bkgs_unc = np.sqrt(tot_bkgs_unc)


                    if left_edge ==  0:
                        custom_binning.append(left_edge)
                        custom_totals.append(curr_integral_signal)
                        done = True
                        break

                    if (tot_integral_signal >= big_bin_counter*bin_content_goal):
                        print(f"Check out the background unc/tot {integral_bkgs_unc/integral_bkgs}")
                        print(f"Check the backgrounds for empty bins")
                        for bkg_name in bkg_names:
                            print(f"{bkg_name} has {integrate_uproot(hists[bkg_name], left_edge, right_edge)[0]} entries")

                        print(f"And the total bkg / tot bkg unc is {tot_bkgs} / {tot_bkgs_unc} = {tot_bkgs/tot_bkgs_unc}")


                    if (tot_integral_signal >= big_bin_counter*bin_content_goal) and (np.max(integral_bkgs_unc/integral_bkgs) <= 0.2) and (tot_bkgs_unc/tot_bkgs <= 0.2):       
                        custom_binning.append(round(hists[signal_name].axis().edges()[left_edge], 2))
                        custom_totals.append(curr_integral_signal)
                        print("At bin ", left_edge, " integral is ", custom_totals[-1])
                        right_edge = left_edge-1
                        big_bin_counter+=1
                        break

            print("Found the bin edge set, it is ", custom_binning)
            print("With totals per bin as ", custom_totals)

            if custom_totals[-1] < 0.2*bin_content_goal:
                print("Final bin was very very small, auto merged!")
                del custom_binning[-2]
                print("New binning is ", custom_binning)

            custom_binning.sort()
            bins_dict[dict_key] = custom_binning

            print("Final dict binnings are ", bins_dict)



            json_dict = convert_to_json(bins_dict, params)

            with open(f'tmp_binning_{sig_string}.json', 'w') as f:
                json.dump(json_dict, f)


            #Create the datacards
            ps_call(f"python3 ../dc_make/create_datacards.py --input /eos/user/d/daebi/bbWW_shape_files --output tmp_datacards_{sig_string} --config ../config/x_hh_bbww_run3.yaml --hist-bins ./tmp_binning_{sig_string}.json --param_values {mass_list}", shell=True)
            #Find out where the limits will be saved
            output = ps_call(f"law run MergeResonantLimits --version dev --datacards 'tmp_datacards_{sig_string}/*.txt' --print-output 0,False", shell=True, catch_stdout=True, split="\n")[1]
            limit_file_name = output[-2]
            #Run the limit calculation
            ps_call(f"law run MergeResonantLimits --version dev --datacards 'tmp_datacards_{sig_string}/*.txt' --remove-output 3,a,y", shell=True)

            #Now load the output numpy array and get value [1] (mass, val, up1, down1, up2, down2)
            print(np.load(limit_file_name).files)
            data = np.load(limit_file_name)
            exp_limit = data[data.files[0]][0][1]

            limit_history[dict_key]['nBins'].append(nbins)
            limit_history[dict_key]['limits'].append(exp_limit)
            limit_list.append(exp_limit)

            #If adding a bin does not improve (limit stays the same) then just move on
            if (exp_limit < best_limit) or (nbins == 1):
                best_limit = exp_limit
                best_bins = bins_dict.copy()

            else:
                print(f"Did not improve limits at nbins {nbins}, best limit was {best_limit}, and current was {exp_limit}")
                #Need to replace the bin dict because it now has the 'worse' new binning
                bins_dict = best_bins.copy()
                break





    print("Finished the limit scan, lets look at the history")
    print(limit_history)

    print(limit_list)
    plt.plot(limit_list)
    #plt.show()
    plt.yscale('log')
    plt.savefig(f"binopt_history_{sig_string}.pdf")
    plt.close()




if __name__ == '__main__':
        file_dir = os.path.dirname(os.path.abspath(__file__))
        pkg_dir = os.path.dirname(file_dir)
        base_dir = os.path.dirname(pkg_dir)
        pkg_dir_name = os.path.split(pkg_dir)[1]
        if base_dir not in sys.path:
            sys.path.append(base_dir)
        __package__ = pkg_dir_name

        from RunKit.run_tools import ps_call

        shapes_dir = "/eos/user/d/daebi/bbWW_shape_files/Run3_2022/"

        masslist = [ 250, 260, 270, 280, 300, 350, 450, 550, 600, 650, 700, 800 ]
        for mass in masslist:
            optimize_binning(shapes_dir+f"shape_m{mass}.root", {'MX': mass}, f'ggRadion_HH_bbWW_M{mass}', mass)
