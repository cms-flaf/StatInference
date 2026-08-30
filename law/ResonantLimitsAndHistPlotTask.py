import law
import luigi

from FLAF.Analysis.tasks import HistPlotTask

from .StatInferenceTask import StatInferenceTask
from .ResonantLimitsTask import ResonantLimitsTask


class ResonantLimitsAndHistPlotTask(StatInferenceTask):
    """The limits, and the plots of the histograms that went into them, from one
    `law run`.

    Depends on ResonantLimitsTask, not PlotResonantLimitsTask: this task is about getting
    the limits together with their input histograms, and the limit plots the configuration
    declares in `limit_plots` are asked for separately, by running PlotResonantLimitsTask.
    """

    workflow = luigi.Parameter(default=law.parameter.NO_STR)

    def get_plottable_eras(self):
        """Real production periods, for FLAF's HistPlotTask.

        Deliberately not the same set as ResonantLimitsTask.get_eras(), which returns the
        eras that carry their own datacards -- for a configuration with an `era_groups:`
        entry that is the meta-era, and this is its members. The two are complements:
        limits are set on the combination, while the input histograms can only be plotted
        for periods that really exist, since HistPlotTask resolves `period` against
        config/<era>/.

        So a group era is replaced by its members rather than dropped. Filtering it out
        and stopping there left a configuration whose `eras:` names only the group -- the
        mainline 2D one does -- with nothing to plot at all, and this task then reported
        success having required no HistPlotTask.
        """
        era_groups = self.get_era_groups()
        eras = self.get_config_data().get("eras", [self.period])
        plottable = []
        for e in eras:
            for real_era in era_groups.get(e, [e]):
                # A configuration may list a group and its members both; each period is
                # plotted once.
                if real_era not in plottable:
                    plottable.append(real_era)
        return plottable

    def requires(self):
        reqs = [ResonantLimitsTask.req(self)]
        for e in self.get_plottable_eras():
            # The histograms to plot are the ones the cards were built from, which live
            # under hists_version -- req() would otherwise copy this task's own version
            # and schedule a fresh production of them under it.
            reqs.append(
                HistPlotTask.req(self, period=e, version=self.input_hists_version)
            )
        return reqs

    def output(self):
        return self.local_target("dummy.txt")

    def run(self):
        self.output().touch()
