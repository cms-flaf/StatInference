import os, sys, ROOT
import yaml

from inference.dhi.scripts.rename_processes import rename_processes
from inference.dhi.scripts.remove_processes import remove_processes

#These lines ensure that the script can import sibling or parent modules/packages, even when executed as a standalone script, by adjusting the Python path and package context dynamically. useful when running scripts from the command line.
file_dir = os.path.dirname(os.path.abspath(__file__))
pkg_dir = os.path.dirname(file_dir)
base_dir = os.path.dirname(pkg_dir)
pkg_dir_name = os.path.split(pkg_dir)[1]
if base_dir not in sys.path:
    sys.path.append(base_dir)
__package__ = pkg_dir_name

input_binning_opt_config = os.path.join(os.environ["ANALYSIS_PATH"], "StatInference", "bin_opt", "bin_optimization.yaml")
with open(input_binning_opt_config, "r") as f:
    input_binning_opt_config_dict = yaml.safe_load(f)

def check_negative_yields(datacard, category, channel):
    for file in os.listdir(input_binning_opt_config_dict["input"]["directory_to_rename_processes"]):
        if (category in file) and (channel.lower() in file or channel in file) and file.endswith(".root"):
            shape_file = ROOT.TFile.Open(os.path.join(input_binning_opt_config_dict["input"]["directory_to_rename_processes"], file))
            print(f"shape file: {shape_file.GetName()}")
    remove_processes_rules = open(f"remove_processes_for_{category}_{channel}_{input_binning_opt_config_dict['input']['era']}.txt", "w")
    remove_processes_list = []

    for key in shape_file.GetListOfKeys():
        hist = key.ReadObj()
        if hist.InheritsFrom("TH1"):
            if hist.Integral() == 0:
                print(f"Zero yield found for histogram {hist.GetName()}")
                remove_processes_rules.writelines(f"{hist.GetName()}\n")
                remove_processes_list.append(hist.GetName())
            elif hist.Integral() < 0:
                print(f"Negative yield found for histogram {hist.GetName()}")
                remove_processes_rules.writelines(f"{hist.GetName()}\n")
                remove_processes_list.append(hist.GetName())
            else:
                continue
        elif hist.InheritsFrom("TDirectoryFile"):
            shape_file.cd(hist.GetName())
            for subkey in ROOT.gDirectory.GetListOfKeys():
                subhist = subkey.ReadObj()
                if subhist.Integral() == 0:
                    print(f"Zero yield found for histogram {subhist.GetName()}")
                    remove_processes_rules.writelines(f"{subhist.GetName()}\n")
                    remove_processes_list.append(subhist.GetName())
                elif subhist.Integral() < 0:
                    print(f"Negative yield found for histogram {subhist.GetName()}")
                    remove_processes_rules.writelines(f"{subhist.GetName()}\n")
                    remove_processes_list.append(subhist.GetName())
                else:
                    continue
            shape_file.cd()
        else:
            continue
    print("remove_processes_list:", remove_processes_list)
    shape_file.Close()
    if len(remove_processes_list) != 0:
        # new_datacard = remove_processes(datacard, remove_processes_list)
        remove_processes(datacard, remove_processes_list)
    else:
        print('No negative or zero yields processes to remove.')
        # new_datacard = datacard
    remove_processes_rules.close()
    # return new_datacard

def copy_processes(category, channel):
    for file in os.listdir(input_binning_opt_config_dict["input"]["directory_to_rename_processes"]):
        if (category in file) and (channel.lower() in file or channel in file) and file.endswith(".root"):
            shape_file = ROOT.TFile.Open(os.path.join(input_binning_opt_config_dict["input"]["directory_to_rename_processes"], file), "UPDATE")
    
    list_of_hists_not_to_copy = []
    for key in shape_file.GetListOfKeys():
        hist = key.ReadObj()

        if hist.InheritsFrom("TH1"):
            list_of_hists_not_to_copy.append(hist.GetName())

    for key in shape_file.GetListOfKeys():
        subdir = key.ReadObj()

        if subdir.InheritsFrom("TDirectoryFile"):
            shape_file.cd(subdir.GetName())
            print(f"In directory {subdir.GetName()} of file {shape_file.GetName()}")
            for subkey in subdir.GetListOfKeys():
                if subkey.ReadObj().InheritsFrom("TH1"):
                    # if subkey.GetName() not in list_of_hists_not_to_copy:
                    if subkey.GetName() == "data_obs":
                        print(f"copying histogram {subkey.GetName()} in shape file {shape_file}")
                        hist_to_copy = ROOT.gDirectory.Get(subkey.GetName())
                        shape_file.cd()
                        hist_to_copy.Write()
                        shape_file.cd(subdir.GetName())
                        # hist_to_copy.Reset()
                    else:
                        continue
                else:
                    continue
    shape_file.Close()

            
if __name__ == "__main__":
    rename_processes_rules = open("my_rules.txt", "r")

    categories = []
    for cat_entry in input_binning_opt_config_dict["input"]["categories_ParameterOfInterest"].split(','):
        category, poi = cat_entry.split(':')
        categories.append(category)
    for channel in input_binning_opt_config_dict["input"]["channels"]:
        for category in categories:
            print("category:", category, "channel:", channel)
            for file in os.listdir(input_binning_opt_config_dict["input"]["directory_to_rename_processes"]):

                if (category in file) and (channel.lower() in file or channel in file) and file.endswith(".txt"):
                    print(f"remove Processing file: {file}")
                    file_path = os.path.join(input_binning_opt_config_dict["input"]["directory_to_rename_processes"], file)
                    check_negative_yields(file_path, category, channel)
                    print(f"rename Processing file: {file}")
                    rename_processes(file_path, rename_processes_rules)
                    print(f"copy Processing file: {file}")
                    copy_processes(category, channel)

    rename_processes_rules.close()