# StatInference

## How to run binning optimisation on lxplus
Open two separate LXPLUS terminals: one for **server side** and the other for **worker side**.
Set up the environment and proxy in analysis area as usual on both terminals.
```sh
source env.sh
voms-proxy-init -voms cms -rfc -valid 192:00
cd StatInference
```
Configure the `bin_opt/bin_optimization.yaml` file for the intended analysis (resonant or nonresonant), channels, era, etc.

Get the server side script running first (it creates worker output directories) and takes about a few seconds. Then move on to the worker side terminal and run the worker side script there.
### Server side
Better to run on screen. If screen does not work, you will need to have this server script below running on an active lxplus terminal.

Example:
For the current version of code, you will need to run the following script separately for each channel.
```sh
python3 bin_opt/optimize_channel.py --channel tauTau
```

### Worker side
Running on screen is not required, jobs are submitted on HTCondor.

Example:
for HTCondor jobs, run the following:
```sh
python3 bin_opt/submitLimitWorkers.py
```
The worker script can also be run locally as following:
```sh
python3 bin_opt/rebinAndRunLimitsWorker.py --channel tauTau
```

If runtime is longer than one day, both server and workers should be restarted.
Workers directory should be removed before restarting server/workers.

E.g.
```sh
rm -r output/binning_v2/tauTau_2018/workers
```
