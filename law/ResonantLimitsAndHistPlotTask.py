import law
import luigi

from FLAF.Analysis.tasks import HistPlotTask

from .StatInferenceTask import StatInferenceTask
from .PlotResonantLimitsTask import PlotResonantLimitsTask


class ResonantLimitsAndHistPlotTask(StatInferenceTask):
    """Everything an analyst wants from one `law run`: the limits, the limit plots, and
    the plots of the histograms that went into them.

    Depends on PlotResonantLimitsTask rather than ResonantLimitsTask -- the plots are the
    product worth asking for, and PlotResonantLimitsTask pulls the limits in behind it, so
    naming the limits here as well would only be a second path to the same task.
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
        """
        era_groups = self.get_era_groups()
        eras = self.get_config_data().get("eras", [self.period])
        return [e for e in eras if e not in era_groups]

    def requires(self):
        reqs = [PlotResonantLimitsTask.req(self)]
        for e in self.get_plottable_eras():
            reqs.append(HistPlotTask.req(self, period=e))
        return reqs

    def output(self):
        return self.local_target("dummy.txt")

    def run(self):
        self.output().touch()
