# StatInference

## Two ways bins get decided

Two independent things in this repository choose the binning of the shapes that go into
the datacards. They do not share code, and an analysis uses one or the other.

**1. In-chain rebinning — `dc_make/hist_rebin_2d.py`.** A production step, run by
`HistRebinTask` in `law/HistRebinTask.py`: it derives the bin edges from the shapes themselves
(no fits, no limits) and writes the rebinned histograms the datacards are built from, so
`CreateDatacardsTask` cannot run without it. Its knobs come from the `binning:` block of
the datacard configuration — see the annotated block in the HH→bbWW Run 3 configuration,
`config/Datacards/x_hh_bbww_DL_run3.yaml` in the analysis area. It cuts each base
category (`SR/res2b`) into per-slice categories whose names come from that block's
`category_pattern`, e.g. `{base_category}_dnn{slice_idx}` -- the pattern is the
analysis's choice, and nothing in this repository assumes the sliced axis is a DNN
score. `common/tools.py:CategoryNaming` both writes those names and parses them back
from the one pattern, and every consumer of the configuration's `categories` list
expands them through it, so the datacard bins and the written histograms cannot fall out
of step.

**2. Offline binning optimisation — `bin_opt/`.** A search harness, documented below: it
builds candidate binnings, runs limits with combine for each, and ranks them. Its product
is a `hist_bins` JSON, applied at datacard time by `dc_make/binner.py`. It is not part of
the datacard chain and exposes no importable API — every module is a script driven by
`bin_opt/bin_optimization.yaml`.

Which one an analysis is on is visible in its configuration: the in-chain path has a
`binning:` block and leaves `hist_bins` unset (as HH→bbWW's `config/global.yaml` does),
while the `bin_opt` path sets `hist_bins` and has no `binning:` block (as
`config/x_hh_bbtautau_run2.yaml` does).

That difference is also what decides how `categories:` is read. With a `binning:` block
(or an explicit `--n-dnn-slices`), `dc_make` expands the listed base categories into the
sliced names above, because those are the directories `HistRebinTask` wrote. Without one,
the categories are used exactly as listed — which is what an input file that is already
1D and binned contains. Pass `--n-dnn-slices 0` to force the latter for a configuration
that does carry a `binning:` block.

In the law chain the same block decides the graph: with no `binning:` block
`HistRebinTask` drops out entirely and `CreateDatacardsTask` reads the merged histograms
straight from the `Hists_merged` tree, since there is no step in between.

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
