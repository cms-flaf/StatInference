#!/bin/bash
# Run the 2D->1D rebinning by hand, the way PreprocessShapesTask runs it.
#
# The chain normally does this for you: the datacard configuration's `preprocess:` block
# names this script, and PreprocessShapesTask supplies --input, --output and --era. This
# is for trying a binning out without running the chain -- point --output somewhere
# scratch and look at what comes out.
set -e

cd "$(dirname "$0")/../.."
source env.sh >/dev/null 2>&1

VERSION=v2605a_DNNOutputs_v3                      # merged 2D histograms to read
CONFIG=config/Datacards/x_hh_bbww_DL_run3.yaml
KNOBS=config/Datacards/binning_2d.yaml
ERA=${1:-Run3_Early}                              # any era of the config's `eras:` list
OUT=${2:-/tmp/$USER/rebin_$ERA}

python3 -u StatInference/bin_opt_2d/rebin_2d.py \
    --input "$ANALYSIS_BIG_DATA_PATH/$VERSION/Hists_merged" \
    --output "$OUT" \
    --config "$CONFIG" \
    --binning-config "$KNOBS" \
    --era "$ERA"

echo "########## DONE ##########"
echo "shapes:  $OUT/<source era>/<variable>/<variable>.root"
echo "binning: $OUT/binning.json"
