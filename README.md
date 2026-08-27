# StatInference

## Two ways bins get decided

Two independent things in this repository choose the binning of the shapes that go into
the datacards. They do not share code, and an analysis uses one or the other. **Neither is
part of the datacard chain**: both are offline pre-steps that produce input the chain then
consumes, so `CreateDatacardsTask` reads whatever shapes it is pointed at and does no
binning of its own.

**1. Shape-driven 2D→1D rebinning — `bin_opt_2d/rebin_2d.py`.** Derives the bin edges from
the shapes themselves (no fits, no limits): it cuts the x axis of each 2D input into slices
that become datacard categories, and rebins y into the bins of each slice. Its knobs come
from `bin_opt_2d/binning.yaml`; `bin_opt_2d/call_rebin_2d.sh` is a worked example.

It writes two things next to each other:

- the rebinned shapes, in the same `<era>/<variable>/<variable>.root` layout FLAF's
  `HistMergerTask` produces — so the chain consumes them by pointing `--hists-version` at
  that production, and nothing downstream is aware a rebinning happened;
- `binning.json`, recording the edges it chose. Passing it back as `--binning` replays
  that binning exactly instead of re-deriving one, which is how a production is frozen and
  reproduced.

Each base category (`SR/res2b`) becomes per-slice categories named by the
`category_pattern` knob, e.g. `{base_category}_dnn{slice_idx}`. The pattern is the
analysis's choice — nothing here assumes the sliced axis is a DNN score.
`common/tools.py:CategoryNaming` both writes those names and parses them back from that one
pattern, so a datacard configuration reading these shapes lists the sliced names in
`categories:` and repeats the same `category_pattern`, which is how the per-category limits
group the slices of one base category back together.

**2. Limit-driven binning optimisation — `bin_opt/`.** A search harness, documented below:
it builds candidate binnings, runs limits with combine for each, and ranks them. Its product
is a `hist_bins` JSON, applied at datacard time by `dc_make/binner.py`. Every module is a
script driven by `bin_opt/bin_optimization.yaml`, and it exposes no importable API.

Which one an analysis is on is visible in its configuration: the `bin_opt` path sets
`hist_bins` (as `config/x_hh_bbtautau_run2.yaml` does), while an analysis reading rebinned
shapes leaves `hist_bins` unset and instead lists the sliced `categories:` and the
`category_pattern` that names them (as HH→bbWW's `config/Datacards/x_hh_bbww_DL_run3.yaml`
does). An analysis using neither simply lists the categories its input already has.

`categories:` is always taken verbatim — it is the set of directories that exist in the
input shapes, and nothing in `dc_make` derives or expands it.

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
