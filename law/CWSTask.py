import law
import luigi
import os
import shlex
import yaml

from FLAF.RunKit.run_tools import ps_call
from FLAF.run_tools.law_customizations import era_csv

from .StatInferenceTask import StatInferenceTask


class CWSTask(StatInferenceTask):
    """Run entries of a CWS analysis config on the datacards of this version, one era at a time.

    CWS (cms-analysis/general/cws, the `cws` submodule of the analysis) is the successor of
    dhi. Its tasks share their names with dhi's, so they live under a law config and task
    index of their own and run on the law CWS pins: they cannot be required from this law
    process, and are started as `cws run` in the environment FLAF's `cwsEnv` alias sets up.

    The analysis config names no era and no path. For each era of `--era` this task exports

    - FLAF_CWS_ERA: the era, which the config maps to its CWS campaign (`flaf.campaigns`);
    - FLAF_CWS_VERSION: `<era>_<version>`, the CWS version, so each era has its own store;
    - FLAF_CWS_DATACARDS: the local mirror of `<fs_StatInference>/<version>/Datacards/<era>`;

    which the config reads back with ``${oc.env:...}``. Each entry of `--cws-task-config` (the
    entries under the config's `tasks:`, all of them by default) is one `cws run`, whose plots
    CWS collects (`--output-dir`) into `<fs_StatInference>/<version>/CWS/<era>/<entry>/`.

    The datacards are inputs: they are read as they are, not built first. fs_StatInference
    falls back to fs_default.
    """

    cws_config = luigi.Parameter(
        default="",
        description="the CWS analysis config, relative to the analysis root or absolute; "
        "defaults to global.yaml's CWS.config",
    )
    cws_task_config = law.CSVParameter(
        default=(),
        description="the entries of the config's `tasks:` to run, comma-separated; all of "
        "them by default",
    )
    cws_workers = luigi.IntParameter(
        default=0,
        significant=False,
        description="luigi workers of each `cws run` (the points of a scan run in parallel); "
        "0 leaves it to the config's `workers`, else 1",
    )
    cws_args = luigi.Parameter(
        default="",
        significant=False,
        description="further options for `cws run`, as one string; write it "
        "--cws-args='--opt value' since it starts with a dash",
    )
    redraw = luigi.BoolParameter(
        default=False,
        significant=False,
        description="discard the existing CWS plots and draw them again "
        "(--remove-output 0,a,y on each `cws run`; no fit is recomputed)",
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Fail before luigi schedules anything: a requested entry the config does not
        # have, or an era it has no campaign for.
        self.cws_entries()
        for era in self.era:
            self.cws_campaign(era)

    # ---- configuration ---------------------------------------------------------------

    def cws_config_path(self):
        path = self.cws_config or self.common(
            lambda setup: (setup.global_params.get("CWS") or {}).get("config"),
            "CWS.config",
        )
        if not path:
            raise ValueError(
                f"{self.get_task_family()}: no CWS analysis config; set CWS.config in "
                "global.yaml or pass --cws-config"
            )
        path = os.path.expandvars(path)
        return path if os.path.isabs(path) else os.path.join(self.ana_path(), path)

    def cws_config_data(self):
        """The analysis config as written. Only the plain keys are read here (`tasks:` and
        `flaf.campaigns`); the ${oc.env:...} interpolations are left to CWS's Hydra."""
        if getattr(self, "_cws_config_data", None) is None:
            with open(self.cws_config_path(), "r") as f:
                self._cws_config_data = yaml.safe_load(f) or {}
        return self._cws_config_data

    def cws_entries(self):
        tasks = list((self.cws_config_data().get("tasks") or {}).keys())
        if not self.cws_task_config:
            if not tasks:
                raise ValueError(f"{self.cws_config_path()} has no `tasks:` to run")
            return tasks
        unknown = [e for e in self.cws_task_config if e not in tasks]
        if unknown:
            raise ValueError(
                f"{self.get_task_family()}: {', '.join(unknown)} not under `tasks:` of "
                f"{self.cws_config_path()}; available: {', '.join(tasks) or 'none'}"
            )
        return list(self.cws_task_config)

    def cws_campaign(self, era):
        campaigns = (self.cws_config_data().get("flaf") or {}).get("campaigns") or {}
        if era not in campaigns:
            raise ValueError(
                f"{self.get_task_family()}(era={era_csv(self.era)}): {era} has no CWS "
                f"campaign; add it to `flaf.campaigns` of {self.cws_config_path()} (a key of "
                "plotting.campaign_labels in cws/cws/conf/general_config.yaml)"
            )
        return str(campaigns[era])

    def cws_version(self, era):
        return f"{era}_{self.version}"

    # ---- law ---------------------------------------------------------------------------

    def datacards_target(self, era):
        return self.output_dir_target(self.version, "Datacards", era)

    def output(self):
        return {
            era: {
                entry: self.output_dir_target(self.version, "CWS", era, entry)
                for entry in self.cws_entries()
            }
            for era in self.era
        }

    def stage_cards(self, era):
        """Mirror an era's datacards to a local directory, returned: combine needs real
        local files. Refreshed on every run, so a rebuilt card is never read stale. Not
        StatInferenceTask.datacards_dir(), which is laid out by StatInference card."""
        source = self.datacards_target(era)
        if not source.exists():
            raise RuntimeError(
                f"{self.get_task_family()}: no datacards for {era} at {source.path}"
            )
        local = law.LocalDirectoryTarget(
            os.path.join(self.ana_data_path(), self.version, "CWS_datacards", era)
        )
        if local.exists():
            local.remove()
        local.touch()
        source.copy_to_local(local)
        return local.abspath

    def cws_env(self, era, datacards_dir):
        """The environment of FLAF's `cwsEnv` alias (CWS's law config, task index and pinned
        law), plus the variables the analysis config reads."""
        env = dict(os.environ)
        cws_base = env.get("CWS_BASE")
        if not cws_base:
            raise RuntimeError(
                "CWS_BASE is not set: CWS is set up by env.sh when HH_CWS_PATH points at the "
                "cws submodule"
            )
        env["LAW_HOME"] = os.path.join(cws_base, ".law")
        env["LAW_CONFIG_FILE"] = os.path.join(cws_base, "law.cfg")
        env["PYTHONPATH"] = ":".join(
            p
            for p in (
                os.path.join(cws_base, "modules", "law"),
                env.get("PYTHONPATH", ""),
            )
            if p
        )
        env["PATH"] = (
            os.path.join(cws_base, "modules", "law", "bin") + ":" + env["PATH"]
        )
        env["FLAF_CWS_ERA"] = era
        env["FLAF_CWS_VERSION"] = self.cws_version(era)
        env["FLAF_CWS_DATACARDS"] = datacards_dir
        return env

    def cws_command(self, entry, output_dir):
        cmd = [
            os.path.join(os.environ["CWS_BASE"], "bin", "cws"),
            "run",
            "--analysis",
            self.cws_config_path(),
            "--task-config",
            entry,
            "--output-dir",
            output_dir,
        ]
        if self.cws_workers > 0:
            cmd += ["--workers", str(self.cws_workers)]
        if self.redraw:
            cmd += ["--remove-output", "0,a,y"]
        cmd += shlex.split(self.cws_args)
        return cmd

    def run(self):
        outputs = self.output()
        for era in self.era:
            env = self.cws_env(era, self.stage_cards(era))
            for entry in self.cws_entries():
                target = outputs[era][entry]
                with target.localize("w") as local_output:
                    # cwd is the analysis root, so a relative path in the config resolves
                    # the way it does for a `cwsEnv cws run` started there.
                    ps_call(
                        self.cws_command(entry, local_output.abspath),
                        env=env,
                        cwd=self.ana_path(),
                        verbose=1,
                    )
                    # CWS only warns when there is nothing to collect, and luigi's default
                    # retcodes exit 0 on a failed task; the collected files are the evidence.
                    if not os.listdir(local_output.abspath):
                        raise RuntimeError(
                            f"CWS entry '{entry}' ({era}) collected nothing into "
                            f"{local_output.abspath}; see the CWS output above for the task "
                            "that failed"
                        )
