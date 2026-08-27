#!/bin/bash
# Rebin the 2D DNN-vs-HME shapes into the sliced 1D shapes the datacards are built from.
#
# A standalone pre-step: the datacard chain does not run this and does not know it
# happened. It writes the same "<era>/<variable>/<variable>.root" layout HistMergerTask
# produces, so the chain consumes the result by pointing --hists-version at $OUT_VERSION.
set -e

cd "$(dirname "$0")/../.."
source env.sh >/dev/null 2>&1

IN_VERSION=v2605a_DNNOutputs        # the merged histograms to read
OUT_VERSION=v2605a_rebinned         # the production this writes
CONFIG=config/Datacards/x_hh_bbww_DL_run3.yaml
BINNING=StatInference/bin_opt_2d/binning.yaml

IN=$ANALYSIS_BIG_DATA_PATH/$IN_VERSION/Hists_merged
OUT=$ANALYSIS_BIG_DATA_PATH/$OUT_VERSION/Hists_merged

# Every sub-era is rebinned with edges discovered from all four summed: the combined fit
# is what the binning is for, and a single era's statistics cannot support the same bins.
ERAS="Run3_2022,Run3_2022EE,Run3_2023,Run3_2023BPix"
for ERA in ${ERAS//,/ }; do
    echo "########## $ERA ##########"
    python3 -u StatInference/bin_opt_2d/rebin_2d.py \
        --input "$IN" \
        --output "$OUT" \
        --config "$CONFIG" \
        --binning-config "$BINNING" \
        --era "$ERA" \
        --discovery-eras "$ERAS"
done

echo "########## DONE ##########"
echo "Rebinned shapes:  $OUT"
echo "Binning recorded: $OUT/binning.json"
echo
echo "Build datacards from them with:"
echo "  law run ResonantLimitsTask --version <v> --period Run3_2022 \\"
echo "      --hists-version $OUT_VERSION"
echo
echo "To reproduce this exact binning later, add:  --binning $OUT/binning.json"
