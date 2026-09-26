import law
import luigi
import shutil
import tempfile

from dhi.tasks.resonant import MergeResonantLimits

from .StatInferenceTask import StatInferenceTask
from .ResonantLimitsTask import (
    ResonantLimitsTask,
    cards_by_mass,
    combine_cards_per_mass,
    publish_dir,
)


class CombinedResonantLimitsTask(StatInferenceTask):
    """The limit of a combination configuration (`members:`): each member's headline
    cards, combined per mass and labelled by member. Same-named nuisances correlate.
    """

    # Not a workflow itself; see ResonantLimitsTask.
    workflow = luigi.Parameter(default=law.parameter.NO_STR)

    def requires(self):
        if not self.is_combination():
            raise RuntimeError(
                f"{self.datacard_config_path()} declares no 'members', so there is "
                "nothing to combine; run ResonantLimitsTask for a single configuration."
            )
        return {
            name: self.member_req(name, ResonantLimitsTask) for name in self.members()
        }

    def store_parts(self):
        return (*self.version_parts(), self.__class__.__name__)

    def output(self):
        return {
            "limits": self.local_target("limits.npz"),
            "datacards": law.LocalDirectoryTarget(self.datacards_dir("combined")),
        }

    @staticmethod
    def check_members(member_cards):
        """The members must share their mass points and no bin."""
        masses = {
            name: sorted(cards_by_mass(cards), key=int)
            for name, cards in member_cards.items()
        }
        reference = next(iter(masses.values()))
        if any(m != reference for m in masses.values()):
            raise RuntimeError(f"the members disagree on the mass points: {masses}")

        owner = {}
        for name, cards in member_cards.items():
            with open(cards[0]) as f:
                # the first "bin" line lists each bin once (the second, per process)
                bins = next(line.split()[1:] for line in f if line.startswith("bin "))
            for b in bins:
                if owner.setdefault(b, name) != name:
                    raise RuntimeError(
                        f"bin '{b}' is in both '{owner[b]}' and '{name}'; the members' "
                        "channels and categories must not overlap"
                    )

    def run(self):
        member_cards = {}
        for name, req in self.requires().items():
            cards = req.headline_cards()
            if not cards:
                raise RuntimeError(f"member '{name}' has no datacards to combine")
            member_cards[name] = cards
        self.check_members(member_cards)

        with tempfile.TemporaryDirectory() as staging:
            combine_cards_per_mass(list(member_cards.items()), staging, self.cmssw_env)
            combined = publish_dir(staging, self.output()["datacards"])

        merged = yield MergeResonantLimits(
            version=self.version, datacards=tuple(sorted(combined))
        )
        self.output()["limits"].parent.touch()
        shutil.copy2(merged.path, self.output()["limits"].path)
