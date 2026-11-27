import glob
import json
import law
import luigi
import os
import shutil
import subprocess
import tempfile

from StatInference.common.tools import importROOT

from .StatInferenceTask import CardTask
from .CreateDatacardsTask import CreateDatacardsTask

# Key of the cross-era combination in the output, next to the era names. Not an era name,
# and deliberately not one a configuration could declare.
COMBINED_KEY = "combined"

# combine's quantileExpected of each expected limit, by the name it is stored under
EXPECTED_QUANTILES = {"2.5": 0.025, "16": 0.16, "50": 0.5, "84": 0.84, "97.5": 0.975}


class NonResonantLimitsTask(CardTask):
    """The upper limits on the signal strength `r` of one non-resonant card (`--card`):
    each of its eras', and the combination of its top-level eras. Without `--card`, the
    limits of every non-resonant card (see CardTask); resonant cards are left out, they
    have ResonantLimitsTask.

    A non-resonant card has no mass to scan, so this is one asymptotic CLs limit per era:
    text2workspace.py with combine's default physics model, whose only POI is `r`, scaling
    every signal process of the card together, and `combine -M AsymptoticLimits` on it.
    No HH model is applied; the limits are on `r` in the units the signal templates are
    normalised to. Blinded (Asimov, a-priori background) unless --unblinded.

    combine runs here rather than through dhi: dhi's UpperLimits always scans a parameter
    other than its POI, and without an HH model there is none to scan.
    """

    # Not a workflow itself. The parameter exists so that --workflow reaches the
    # CreateDatacardsTask below -- law's req() only forwards parameters both tasks declare
    # -- and NO_STR leaves it at its own default. The fits themselves run here, locally.
    workflow = luigi.Parameter(default=law.parameter.NO_STR)

    unblinded = luigi.BoolParameter(
        default=False,
        description="also compute the observed limit; default: blinded, expected only",
    )

    check_card_eras = True

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.card and not self.is_nonresonant(self.card):
            raise ValueError(
                f"{self.get_task_family()}: card '{self.card}' has model parameters, so it "
                "is resonant: its limits are ResonantLimitsTask's"
            )

    def card_has_work(self, card):
        return self.is_nonresonant(card)

    def store_parts(self):
        return super().store_parts() + ("combined",)

    def requires(self):
        if self.all_cards:
            return self.card_tasks()
        return CreateDatacardsTask.req(self, branches=self.card_branches())

    def card_branches(self):
        """The CreateDatacardsTask branches of this card, as a `branches` selection: only
        they are built, so one card's limits neither wait for nor fail with another card.
        () -- the whole workflow -- when this card is every branch, as for a single card.
        """
        workflow = CreateDatacardsTask.req(self, branches=())
        branch_map = workflow.get_branch_map()
        indices = tuple(
            index
            for index, data in branch_map.items()
            if data["card"] == self.this_card
        )
        if not indices:
            raise ValueError(
                f"{self.get_task_family()}: card '{self.this_card}' has no datacards to build"
            )
        return () if len(indices) == len(branch_map) else indices

    def output(self):
        if self.all_cards:
            return {card: task.output() for card, task in self.card_tasks().items()}
        return {
            # {era or "combined": {"expected": {quantile: limit}, "observed": limit or None}}
            "limits": self.local_target("limits.json"),
            "datacards": law.LocalDirectoryTarget(self.datacards_dir(COMBINED_KEY)),
        }

    def stage_datacards(self, era, remote_target):
        """Mirror an era's datacards from fs_StatInference to datacards_dir(era), returned:
        combine needs real local files. Re-staged on every run, so the mirror never goes
        stale when CreateDatacardsTask reruns within the same version."""
        local_dir = self.datacards_dir(era)
        local_target = law.LocalDirectoryTarget(local_dir)
        if local_target.exists():
            local_target.remove()
        part = law.LocalDirectoryTarget(is_tmp=True)
        part.touch()
        try:
            remote_target.copy_to_local(part)
            shutil.copytree(part.abspath, local_dir)
        finally:
            part.remove(silent=True)
        return local_dir

    def era_card(self, era, card_dir, staging):
        """The one datacard of an era: datacard_combined_signals.txt when the card keeps
        every signal in one datacard, else the per-bin datacards of the top-level
        directory combined into one (in `staging`)."""
        cards = sorted(glob.glob(os.path.join(card_dir, "*.txt")))
        if not cards:
            raise RuntimeError(f"no datacard in {card_dir} for era {era}")
        if len(cards) == 1:
            return cards[0]
        return self.combine_cards(
            [(f"b{i}", c) for i, c in enumerate(cards)],
            os.path.join(staging, f"datacard_{era}.txt"),
        )

    def combine_cards(self, named_cards, out_file):
        """combineCards.py of [(name, card)] into `out_file`. The cards are given by
        absolute path, so the shape files they point to stay where they are."""
        cmd = ["combineCards.py"] + [
            f"{name}={os.path.abspath(card)}" for name, card in named_cards
        ]
        with open(out_file, "w") as f:
            subprocess.run(cmd, env=self.cmssw_env, stdout=f, check=True)
        return out_file

    def asymptotic_limits(self, card, label):
        """{"expected": {quantile: limit}, "observed": limit or None} of a datacard."""
        with tempfile.TemporaryDirectory() as work:
            workspace = os.path.join(work, "workspace.root")
            cmd = ["text2workspace.py", os.path.abspath(card), "-o", workspace]
            # A bin with no background makes text2workspace refuse the card. A card opts
            # in to accepting it -- the integration CI, whose few events can leave a
            # channel x category of an era without any background -- by setting
            # allow_no_background; the analysis cards keep the check.
            if self.get_config_data(self.this_card).get("allow_no_background", False):
                cmd.append("--X-allow-no-background")
            subprocess.run(
                cmd,
                env=self.cmssw_env,
                cwd=work,
                check=True,
            )
            cmd = ["combine", "-M", "AsymptoticLimits", workspace, "-n", "_limits"]
            if not self.unblinded:
                cmd += ["--run", "blind"]
            subprocess.run(cmd, env=self.cmssw_env, cwd=work, check=True)
            (result,) = glob.glob(
                os.path.join(work, "higgsCombine_limits.AsymptoticLimits.*.root")
            )
            limits = self.read_limits(result, observed=self.unblinded)
        print(f"{self.get_task_family()}: {label}: {limits}")
        return limits

    @staticmethod
    def read_limits(path, observed=False):
        """The limits combine wrote to `path`. A quantile combine did not compute is None,
        with a warning, as dhi's UpperLimits.load_limits leaves it NaN: combine can exit 0
        after a fit it could not finish -- likely at --test statistics, where the CI only
        has to run the chain through. Raises only when there is no limit tree at all."""
        ROOT = importROOT()
        f = ROOT.TFile.Open(path)
        try:
            tree = f.Get("limit")
            if not tree:
                raise RuntimeError(f"{path} holds no limit tree: the fit failed")
            values = [(round(entry.quantileExpected, 3), entry.limit) for entry in tree]
        finally:
            f.Close()
        expected = {}
        for name, quantile in EXPECTED_QUANTILES.items():
            matches = [v for q, v in values if abs(q - quantile) < 1e-3]
            expected[name] = matches[0] if matches else None
        observed_values = [v for q, v in values if q == -1]
        missing = [f"expected {name}%" for name, v in expected.items() if v is None]
        if observed and not observed_values:
            missing.append("observed")
        if missing:
            print(
                f"WARNING: {path}: combine did not compute the {', '.join(missing)} "
                "limit(s); recorded as null (see combine's output above)"
            )
        return {
            "expected": expected,
            "observed": observed_values[0] if observed else None,
        }

    def run(self):
        if self.all_cards:
            # every card's limits are its requirements; nothing is left to do here
            return
        era_outputs = {
            task.branch_data["era"]: task.output()
            for task in self.requires().get_branch_tasks().values()
            if task.branch_data["card"] == self.this_card
        }
        limits = {}
        with tempfile.TemporaryDirectory() as staging:
            era_cards = {
                era: self.era_card(era, self.stage_datacards(era, output), staging)
                for era, output in era_outputs.items()
            }
            for era, card in era_cards.items():
                limits[era] = self.asymptotic_limits(card, era)

            # The cross-era combination is the one place a group era and its members
            # cannot both appear: the group already *is* their combination, so a card
            # built from both would count those events twice.
            top_level = [e for e in self.get_top_level_eras() if e in era_cards]
            if not top_level:
                raise RuntimeError(
                    "No datacards were found for any top-level era "
                    f"({', '.join(self.get_top_level_eras())}), so there is no combined "
                    "limit to write."
                )
            out_dc_dir = self.output()["datacards"]
            if out_dc_dir.exists():
                out_dc_dir.remove()
            out_dc_dir.touch()
            combined_card = self.combine_cards(
                [(e, era_cards[e]) for e in top_level],
                os.path.join(out_dc_dir.path, "datacard_combined.txt"),
            )
            # With a single top-level era its own limit already is the combination.
            limits[COMBINED_KEY] = (
                limits[top_level[0]]
                if len(top_level) == 1
                else self.asymptotic_limits(combined_card, COMBINED_KEY)
            )

        self.output()["limits"].parent.touch()
        with open(self.output()["limits"].path, "w") as f:
            json.dump(
                {"poi": "r", "blinded": not self.unblinded, "limits": limits},
                f,
                indent=2,
            )
