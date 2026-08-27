import glob
import law
import luigi
import os
import re
import shutil
import subprocess

from dhi.tasks.resonant import MergeResonantLimits

from .StatInferenceTask import StatInferenceTask
from .CreateDatacardsTask import CreateDatacardsTask


class ResonantLimitsTask(StatInferenceTask):
    workflow = luigi.Parameter(default=law.parameter.NO_STR)

    def store_parts(self):
        return (self.version, self.__class__.__name__, "combined")

    def get_eras(self):
        return self.get_top_level_eras()

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

        self.output()["limits"].parent.touch()
        shutil.copy2(limits.path, self.output()["limits"].path)

        out_dc_dir = self.output()["datacards"]
        out_dc_dir.touch()

        masses = set()
        for e, cards in era_cards.items():
            for c in cards:
                m = re.search(r"_(\d+)\.txt$", c)
                if m:
                    masses.add(m.group(1))

        for mass in masses:
            combine_args = []
            for e in eras:
                for c in era_cards[e]:
                    if c.endswith(f"_{mass}.txt"):
                        combine_args.append(f"{e}={c}")
                        break

            if combine_args:
                cmd = ["combineCards.py"] + combine_args
                out_file = os.path.join(out_dc_dir.path, f"combined_{mass}.txt")
                with open(out_file, "w") as f:
                    subprocess.run(cmd, env=self.cmssw_env, stdout=f, check=True)
