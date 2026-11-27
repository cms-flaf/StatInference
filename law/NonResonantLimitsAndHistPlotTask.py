import law
import luigi

from FLAF.Analysis.tasks import HistPlotTask

from .StatInferenceTask import StatInferenceTask
from .NonResonantLimitsTask import NonResonantLimitsTask


class NonResonantLimitsAndHistPlotTask(StatInferenceTask):
    """The limits of the non-resonant cards, and the plots of the histograms that went into
    them, from one `law run`: the non-resonant counterpart of ResonantLimitsAndHistPlotTask,
    and the task the integration pipeline runs for a non-resonant analysis.

    Its output is a marker written once everything below it ran, so a version it finished
    is not looked at again: a later run of that version over more eras, or with more cards,
    builds and plots nothing new.
    Run HistPlotTask or NonResonantLimitsTask directly for that, or remove the marker first
    (`--remove-output 0,a,y`).
    """

    # Not a workflow itself. The parameter exists so that --workflow reaches the tasks
    # below that are -- law's req() only forwards parameters both tasks declare -- and
    # NO_STR leaves each of those at its own default.
    workflow = luigi.Parameter(default=law.parameter.NO_STR)

    check_card_eras = True

    def requires(self):
        # The histograms to plot are the ones the cards were built from, which live under
        # hists_version -- req() would otherwise copy this task's own version and schedule
        # a fresh production of them under it. One HistPlotTask over every era of --era,
        # which holds every real era behind the cards (and may hold more, plotted too).
        return [
            NonResonantLimitsTask.req(self),
            HistPlotTask.req(self, version=self.input_hists_version),
        ]

    def output(self):
        return self.local_target("dummy.txt")

    def run(self):
        self.output().touch()
