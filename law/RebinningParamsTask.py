import hashlib
import luigi
import os

from StatInference.common.tools import CategoryNaming

from .StatInferenceTask import StatInferenceTask

# Binning knobs whose values change the rebinned shapes. Kept in one place because they
# are propagated down the chain, read from the datacard configuration's `binning:` block,
# and hashed into the output path -- see RebinningParamsTask.binning_parts().
BINNING_PARAMS = (
    "n_slices",
    "category_pattern",
    "slice_var",
    "max_bins_per_slice",
    "min_slice_bkg_sum",
    "min_bin_bkg_each",
    "min_signal",
    "min_slice_bkg_neff",
    "min_bkg_frac",
    "min_bin_bkg_neff",
    "bkg_per_bin",
    "significance_mode",
    "min_slice_bkg_each",
    "min_slice_bkg_each_neff",
)

# Fallbacks for a configuration that declares no `binning:` block. The values that matter
# belong to an analysis, not to this file: what a slice needs to be worth keeping depends
# on the sample sizes and the selection, and a test configuration wants them at zero. Each
# datacard configuration should state its own -- see the annotated block in
# config/Datacards/x_hh_bbww_DL_run3.yaml for what the production numbers are and why.
BINNING_DEFAULTS = {
    "n_slices": 4,
    # Neutral by default: the sliced axis is whatever the analysis put on the 2D x axis,
    # and an analysis that wants its own name for it says so in its `binning:` block.
    "category_pattern": CategoryNaming.default_pattern,
    "slice_var": "x",
    "max_bins_per_slice": 10,
    "min_slice_bkg_sum": 1.0,
    "min_bin_bkg_each": 0.01,
    "min_signal": 0.5,
    "min_slice_bkg_neff": 4.0,
    "min_bkg_frac": 0.05,
    "min_bin_bkg_neff": 4.0,
    "bkg_per_bin": 5.0,
    "significance_mode": "asimov",
    "min_slice_bkg_each": 0.01,
    "min_slice_bkg_each_neff": 0.0,
}


class RebinningParamsTask(StatInferenceTask):
    """The binning parameters, declared once for the whole chain.

    They live here rather than on HistRebinTask alone for two reasons. First, `.req()`
    only copies parameters the requesting task also declares, so with them on HistRebinTask
    only there was no way to steer the binning from ResonantLimitsTask -- the top of the
    chain always got the defaults. Second, every task's output path has to agree on which
    binning it describes, which needs the values everywhere.
    """

    # All default to None, meaning "take it from the datacard configuration". Optional*
    # rather than plain Int/Float/Parameter is required, not cosmetic: luigi serialises a
    # plain FloatParameter's None to the string "None" and then raises on float("None")
    # whenever it round-trips the task through to_str_params()/from_str_params().
    n_slices = luigi.OptionalIntParameter(default=None)
    category_pattern = luigi.OptionalParameter(default=None)
    slice_var = luigi.OptionalParameter(default=None)
    max_bins_per_slice = luigi.OptionalIntParameter(default=None)
    min_slice_bkg_sum = luigi.OptionalFloatParameter(default=None)
    min_bin_bkg_each = luigi.OptionalFloatParameter(default=None)
    min_signal = luigi.OptionalFloatParameter(default=None)
    min_slice_bkg_neff = luigi.OptionalFloatParameter(default=None)
    min_bkg_frac = luigi.OptionalFloatParameter(default=None)
    min_bin_bkg_neff = luigi.OptionalFloatParameter(default=None)
    bkg_per_bin = luigi.OptionalFloatParameter(default=None)
    significance_mode = luigi.OptionalParameter(default=None)
    min_slice_bkg_each = luigi.OptionalFloatParameter(default=None)
    min_slice_bkg_each_neff = luigi.OptionalFloatParameter(default=None)

    def rebinning_enabled(self):
        """Whether the configuration asks for the in-chain 2D->1D rebinning step.

        A configuration with no ``binning:`` block describes input that is already 1D and
        binned: HistRebinTask has nothing to do with it, and CreateDatacardsTask reads the
        merged histograms directly. This is the same signal DatacardMaker uses to decide
        whether the categories are sliced, so the two cannot disagree.
        """
        return bool(self.get_config_data().get("binning"))

    def config_binning(self):
        """The datacard configuration's ``binning:`` block, filled in from
        BINNING_DEFAULTS for anything it does not declare."""
        declared = self.get_config_data().get("binning") or {}
        unknown = set(declared) - set(BINNING_PARAMS)
        if unknown:
            raise RuntimeError(
                f"{self.datacard_config_path()}: unknown key(s) in the 'binning' block: "
                f"{sorted(unknown)}. Known keys: {sorted(BINNING_PARAMS)}."
            )
        return {p: declared.get(p, BINNING_DEFAULTS[p]) for p in BINNING_PARAMS}

    def binning_params(self):
        """Effective binning: an explicitly given parameter wins, otherwise the
        configuration's value."""
        config = self.config_binning()
        return {
            p: (getattr(self, p) if getattr(self, p) is not None else config[p])
            for p in BINNING_PARAMS
        }

    def binning_diffs(self):
        """Binning parameters overridden away from what the configuration asks for."""
        config = self.config_binning()
        return {p: v for p, v in self.binning_params().items() if v != config[p]}

    def binning_parts(self):
        """Extra path level identifying a non-default binning, as a tuple to splat.

        law decides completeness from output existence alone, so without this a changed
        binning parameter silently reuses the previous shapes -- the task reports complete
        and the limits are computed from the old edges.

        Empty for the default binning, so nominal output paths (and the
        --multi-datacards globs that point at them) are exactly where they always were,
        and only scan variants gain a directory level. HistRebinTask writes
        binning_params.json into its output so a hash directory can be decoded later.
        """
        diffs = self.binning_diffs()
        if not diffs:
            return ()
        payload = ",".join(f"{k}={diffs[k]!r}" for k in sorted(diffs))
        return ("bin_" + hashlib.sha1(payload.encode()).hexdigest()[:10],)

    def store_parts(self):
        return (
            self.version,
            self.__class__.__name__,
            *self.binning_parts(),
            self.period,
        )

    def datacards_dir(self, era):
        """Local directory holding an era's datacards.

        The cards are produced on fs_default but must be real local files by the time
        combine sees them; ResonantLimitsTask mirrors them here and everything downstream
        (dhi's --multi-datacards globbing, the overlay plots) resolves against this path.
        """
        return os.path.join(
            self.ana_data_path(),
            self.version,
            "Datacards",
            *self.binning_parts(),
            era,
        )
