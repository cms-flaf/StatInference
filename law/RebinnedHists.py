import law
import luigi

from .StatInferenceTask import StatInferenceTask


class RebinnedHists(StatInferenceTask, law.ExternalTask):
    """One rebinned histogram file, as bin_opt_2d/rebin_2d.py wrote it.

    External for the same reason MergedHists is: the rebinning is a standalone pre-step,
    so a missing input is a missing dependency rather than something this chain can make.

    The path carries two eras, and they are not the same thing. `target_era` is the era
    whose binning this is -- the era of the datacard configuration's `eras:` list that
    rebin_2d.py was asked to produce. `period` is the era whose events these are. For a
    standalone era the two match; for a group era every member appears under the group's
    binning, which is why the target era is a directory level rather than part of the
    version: the same period appears under several binnings and they must not collide.
    """

    variable = luigi.Parameter(description="rebinned variable (directory) name")
    target_era = luigi.Parameter(description="era whose binning this file was cut with")

    def output(self):
        return self.remote_target(
            self.input_hists_version,
            "Hists_rebinned",
            self.target_era,
            self.period,
            self.variable,
            f"{self.variable}.root",
            fs=self.fs_HistTuple,
        )
