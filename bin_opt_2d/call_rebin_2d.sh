#!/bin/bash
# Derive the binning for every era of the datacard configuration's `eras:` list.
#
# A standalone pre-step: the datacard chain does not run this. Its product is a
# binning.json per era, which CreateDatacardsTask applies to the merged 2D shapes at
# datacard time via the configuration's `hist_bins`. No intermediate shapes are written,
# so nothing can go stale against the binning that describes it.
set -e

cd "$(dirname "$0")/../.."
source env.sh >/dev/null 2>&1

VERSION=v2605a_DNNOutputs_v3                     # the merged 2D histograms to read
CONFIG=config/Datacards/x_hh_bbww_DL_run3.yaml
BINNING=StatInference/bin_opt_2d/binning.yaml
OUT=config/Datacards/binnings                    # committed alongside the card

IN=$ANALYSIS_BIG_DATA_PATH/$VERSION/Hists_merged

# Every era the configuration lists, including the group. A plain era is binned on its own
# statistics for a standalone limit; a group era on its members' summed statistics. They
# are different binnings and each gets its own file.
ERAS=$(python3 -c "import yaml;print(' '.join(yaml.safe_load(open('$CONFIG'))['eras']))")

for ERA in $ERAS; do
    echo "########## $ERA ##########"
    python3 -u StatInference/bin_opt_2d/rebin_2d.py \
        --input "$IN" --output "$OUT" \
        --config "$CONFIG" --binning-config "$BINNING" --era "$ERA"
done

echo "########## DONE ##########"
echo "Point the configuration's hist_bins at the era you are building cards for, e.g."
echo "  hist_bins: $OUT/Run3_Early/binning.json"
