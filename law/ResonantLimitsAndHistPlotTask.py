import law
import luigi

from FLAF.Analysis.tasks import HistPlotTask

from .RebinningParamsTask import RebinningParamsTask
from .ResonantLimitsTask import ResonantLimitsTask


class ResonantLimitsAndHistPlotTask(RebinningParamsTask):
    workflow = luigi.Parameter(default=law.parameter.NO_STR)

    def get_eras(self):
        """Real, plottable periods only -- meta-era names are excluded since
        HistPlotTask requires an actual production period."""
        era_groups = self.get_era_groups()
        eras = self.get_config_data().get("eras", [self.period])
        return [e for e in eras if e not in era_groups]

    def requires(self):
        reqs = [ResonantLimitsTask.req(self)]
        for e in self.get_eras():
            reqs.append(HistPlotTask.req(self, period=e))
        return reqs

    def output(self):
        return self.local_target("dummy.txt")

    def run(self):
        self.output().touch()
