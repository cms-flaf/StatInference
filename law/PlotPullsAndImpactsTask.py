import glob
import json
import law
import luigi
import os
import re

from string import Template

from dhi.tasks.pulls_impacts import MergePullsAndImpacts, PlotPullsAndImpacts

from .DhiPlotMixin import DhiPlotMixin
from .StatInferenceTask import StatInferenceTask
from .ResonantLimitsTask import ResonantLimitsTask


class PlotPullsAndImpactsTask(DhiPlotMixin, StatInferenceTask):
    """The nuisance pull and impact plots, declared in the datacard configuration.

    Each entry of the config's ``impact_plots`` block becomes one dhi PlotPullsAndImpacts
    task per mass point it lists. Like the limit plots, the results are copied to
    fs_default beside the rest of the chain's products; dhi keeps its own copy under
    inference/data/store, which is what makes a redraw cheap.

    The input is the combined card ResonantLimitsTask already writes -- one per mass,
    every era in it, at ``<data>/<version>/Datacards/combined/combined_<mass>.txt``. Their
    shape paths are absolute, which is what makes them usable here at all: PullsAndImpacts
    runs combine in a temporary directory.

    This is a task of its own rather than something the limit chain drags in because it is
    far more expensive than everything before it. PullsAndImpacts is a per-parameter
    workflow: roughly two combine fits per nuisance per mass point, so a card with 76
    nuisances costs ~150 fits for a single mass.
    """

    workflow = luigi.Parameter(default=law.parameter.NO_STR)

    def requires(self):
        # ResonantLimitsTask is what writes the combined cards this reads.
        return [ResonantLimitsTask.req(self)]

    def output(self):
        return self.output_dir_target(self.version, "ImpactPlots")

    def entry_eras(self, entry):
        """The eras an entry is drawn for.

        An entry with a `glob:` addresses cards inside one era's directory, so it is
        drawn once per era the configuration lists. The default combined card already
        contains every era, so it is drawn once and labelled "combined".
        """
        return self.get_all_eras() if entry.get("glob") else ["combined"]

    def card_for(self, entry, era, mass):
        """The single datacard an entry's plot is built from.

        PullsAndImpacts fits one card, not a glob -- every parameter gets its own
        combine job against that one workspace -- so an entry's `glob:` and mass have to
        resolve to exactly one file.
        """
        name = entry.get("name") or "impacts"
        pattern = entry.get("glob")
        if not pattern:
            # The combined card ResonantLimitsTask writes: one per mass, every era in it.
            base = self.datacards_dir("combined")
            candidates = glob.glob(os.path.join(base, "combined_*.txt"))
        else:
            base = self.datacards_dir(era)
            candidates = glob.glob(
                os.path.join(base, Template(pattern).safe_substitute(ERA=era))
            )

        # Same convention dhi's datacard_pattern uses: the mass is the trailing number.
        matched = [
            p
            for p in candidates
            if (m := re.search(r"_(\d+)\.txt$", os.path.basename(p)))
            and int(m.group(1)) == int(mass)
        ]
        if len(matched) != 1:
            found = sorted(os.path.basename(p) for p in candidates)
            raise RuntimeError(
                f"impact_plots entry '{name}' ({era}, mass {mass}): expected exactly one "
                f"datacard, found {len(matched)}. Looked in {base} for "
                f"{pattern or 'combined_*.txt'}; it holds: {', '.join(found) or 'nothing'}."
            )
        return matched[0]

    def build_plot_spec(self, entry, era, mass):
        """PlotPullsAndImpacts parameters for an entry and mass point.

        Kept as a plain dict so the same values can both construct the task (to learn
        where it will write) and be rendered into a command line (to make it write).
        """
        name = entry.get("name") or "impacts"
        plot_params = dict(entry.get("plot_params") or {})
        self.validate_plot_params(
            PlotPullsAndImpacts, plot_params, f"impact_plots entry '{name}'"
        )

        # The combined cards carry a few hundred autoMCStats bins, so mc_stats puts ~460
        # parameters on the plot. parameters_per_page defaults to -1, meaning a single
        # page, and the result is an unreadable hairline strip rather than an error.
        if (
            plot_params.get("mc_stats")
            and plot_params.get("parameters_per_page", -1) < 1
        ):
            raise RuntimeError(
                f"impact_plots entry '{name}': mc_stats puts every autoMCStats bin on the "
                "plot, and parameters_per_page defaults to -1 (one page), which is "
                "unreadable. Set parameters_per_page (25 works) or drop mc_stats."
            )

        return dict(
            version=self.version,
            datacards=(self.card_for(entry, era, mass),),
            mass=float(mass),
            # PullsAndImpacts is a plain POITask, so hh_model defaults to the
            # non-resonant model_default -- which would fit r together with kl, kt, CV and
            # C2V, and build its own workspace under a hh_model__model_default/ store path
            # rather than reusing the one the resonant chain already made. dhi's own
            # resonant tasks pin hh_model to NO_STR for the same reason
            # (dhi/tasks/resonant.py:35); allow_empty_hh_model is already True on the
            # PullsAndImpacts base, so this is simply saying the analysis is resonant.
            # An entry may still override it through plot_params.
            **{"hh_model": law.NO_STR, **plot_params},
        )

    def report_dropped_parameters(self, entry, era, mass, spec):
        """Warn about nuisances robustHesse removed from the fit.

        robustHesse drops a parameter it cannot invert, logs "Dropping <name> from the
        hessian" and carries on successfully; the nuisance is then simply absent from the
        plot and the merged JSON, with nothing on the plot marking its absence. So the
        card's own nuisance list is the reference: anything in it that the merged JSON
        does not carry was dropped.
        """
        merged = MergePullsAndImpacts(
            **{
                k: v
                for k, v in spec.items()
                if k in self.known_plot_params(MergePullsAndImpacts)
            }
        ).output()
        if not os.path.exists(merged.path):
            return
        with open(merged.path) as f:
            fitted = {p["name"] for p in json.load(f).get("params", [])}

        card_params = set()
        with open(self.card_for(entry, era, mass)) as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2 and parts[1] in ("shape", "lnN"):
                    card_params.add(parts[0])

        dropped = sorted(card_params - fitted)
        if dropped:
            name = entry.get("name") or "impacts"
            print(
                f"WARNING: impact_plots entry '{name}' ({era}, mass {mass}): "
                f"{len(dropped)} "
                f"nuisance(s) are in the card but not in the fit, so they are missing "
                f"from the plot without being marked: {', '.join(dropped)}. "
                "With --method robust this is robustHesse dropping what it could not "
                "invert; the ranking is not complete."
            )

    def run(self):
        entries = self.get_config_data().get("impact_plots", [])
        if not entries:
            raise RuntimeError(
                f"{self.datacard_config_path()} declares no 'impact_plots' block, so "
                "there is nothing to plot."
            )

        with self.output().localize("w") as local_output:
            produced = []
            for entry in entries:
                name = entry.get("name") or "impacts"
                masses = entry.get("masses")
                if not masses:
                    raise RuntimeError(
                        f"impact_plots entry '{name}': no 'masses'. Each one costs about "
                        "two combine fits per nuisance, so they are listed explicitly "
                        "rather than defaulting to every mass in the model."
                    )

                for era in self.entry_eras(entry):
                    for mass in masses:
                        rel = os.path.join(name, era, str(mass))
                        dest_dir = os.path.join(local_output.abspath, rel)
                        os.makedirs(dest_dir, exist_ok=True)

                        spec = self.build_plot_spec(entry, era, mass)
                        basenames = self.draw_plot(
                            PlotPullsAndImpacts,
                            spec,
                            name,
                            f"{era} mass {mass}",
                            dest_dir,
                            extra_args=entry.get("dhi_args") or (),
                        )
                        self.report_dropped_parameters(entry, era, mass, spec)
                        produced.extend(os.path.join(rel, b) for b in basenames)

            with open(os.path.join(local_output.abspath, "impacts.json"), "w") as f:
                json.dump(sorted(produced), f, indent=2)
