import glob
import law
import luigi
import os
import re
import shutil
import subprocess
import tempfile

from dhi.tasks.resonant import MergeResonantLimits

from .StatInferenceTask import StatInferenceTask
from .CreateDatacardsTask import CreateDatacardsTask


class ResonantLimitsTask(StatInferenceTask):
    # Not a workflow itself. The parameter exists so that --workflow reaches the tasks
    # below that are -- law's req() only forwards parameters both tasks declare -- and
    # NO_STR leaves each of those at its own default.
    workflow = luigi.Parameter(default=law.parameter.NO_STR)

    def store_parts(self):
        return (self.version, self.__class__.__name__, "combined")

    def get_eras(self):
        """Every era the configuration lists: each gets its own cards and its own limit."""
        return self.get_all_eras()

    def _create_datacards_req(self, era, **kwargs):
        era_groups = self.get_era_groups()
        if era in era_groups:
            return CreateDatacardsTask.req(
                self, period=era_groups[era][0], meta_era=era, **kwargs
            )
        return CreateDatacardsTask.req(self, period=era, **kwargs)

    def requires(self):
        return [self._create_datacards_req(e, branches=()) for e in self.get_eras()]

    def output(self):
        return {
            "limits": self.local_target("limits.npz"),
            "datacards": law.LocalDirectoryTarget(self.datacards_dir("combined")),
        }

    def stage_datacards(self, era, remote_target):
        """Mirror an era's datacards from fs_default to a stable local path, returned.

        CreateDatacardsTask writes to fs_default, but the cards must be real local files
        by the time combine sees them: MergeResonantLimits shells out to combine, and
        PlotMultipleResonantLimits resolves --multi-datacards by globbing the filesystem
        outside law entirely. A localize() scratch dir will not do -- run() yields, which
        suspends and re-enters it, so the scratch dir is gone before dhi reads anything.

        The mirror path is the one CreateDatacardsTask used to write to directly, so
        existing --multi-datacards invocations keep working unchanged.
        """
        local_dir = self.datacards_dir(era)
        local_target = law.LocalDirectoryTarget(local_dir)
        # Re-staged on every run, including the re-entry after the yield below. The
        # copy is small (~30 MB) and unconditional refresh is what keeps the mirror
        # from going stale when CreateDatacardsTask reruns within the same version.
        if local_target.exists():
            local_target.remove()
        local_target.touch()
        remote_target.copy_to_local(local_target)
        return local_dir

    def run(self):
        datacards = []
        eras = self.get_eras()
        era_cards = {}

        for e in eras:
            create_dc_br0 = self._create_datacards_req(e, branch=0, branches=())
            output_dir = self.stage_datacards(e, create_dc_br0.output())
            cards = glob.glob(os.path.join(output_dir, "*.txt"))
            era_cards[e] = cards
            datacards.extend(cards)

        limits = yield MergeResonantLimits(
            version=self.version, datacards=tuple(datacards)
        )
        print(f"Merged limits: {limits}")

        masses = set()
        for e, cards in era_cards.items():
            for c in cards:
                m = re.search(r"_(\d+)\.txt$", c)
                if m:
                    masses.add(m.group(1))

        # Built into a staging directory and published only once every card is there.
        # Creating the outputs first and filling them as we go left a failed
        # combineCards halfway through looking like a complete task -- law's default
        # completeness is output existence -- so the next run skipped it and everything
        # downstream read a partial set. Replacing the directory rather than writing into
        # it also clears cards for masses no longer in the configuration.
        with tempfile.TemporaryDirectory() as staging:
            # The cross-era combination is the one place a group era and its members
            # cannot both appear: the group already *is* their combination, so a card
            # built from both would count those events twice.
            for mass in sorted(masses):
                combine_args = []
                for e in self.get_top_level_eras():
                    for c in era_cards.get(e, []):
                        if c.endswith(f"_{mass}.txt"):
                            combine_args.append(f"{e}={c}")
                            break

                if combine_args:
                    cmd = ["combineCards.py"] + combine_args
                    out_file = os.path.join(staging, f"combined_{mass}.txt")
                    with open(out_file, "w") as f:
                        subprocess.run(cmd, env=self.cmssw_env, stdout=f, check=True)

            out_dc_dir = self.output()["datacards"]
            if out_dc_dir.exists():
                out_dc_dir.remove()
            out_dc_dir.touch()
            for name in sorted(os.listdir(staging)):
                shutil.copy2(
                    os.path.join(staging, name), os.path.join(out_dc_dir.path, name)
                )

        self.output()["limits"].parent.touch()
        shutil.copy2(limits.path, self.output()["limits"].path)
