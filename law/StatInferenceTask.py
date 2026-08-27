import contextlib
import law
import luigi
import os
import shutil
import yaml

from string import Template
from FLAF.run_tools.law_customizations import Task


class StatInferenceTask(Task):
    """Shared datacard-configuration access for the limit-setting chain.

    Everything downstream of the merged histograms -- rebinning, datacards, limits, limit
    plots -- is driven by the datacard configuration named in global.yaml's
    ``StatInference.config``, not by global.yaml's own variable lists. This base class is
    the single place that file is read, so the four tasks below cannot disagree about
    which eras, masses or variables the analysis consists of.
    """

    # The Hists_merged tree to read. Separate from `version`, which names what this chain
    # *writes*: re-binning or re-fitting under a new version must not require the input
    # histograms to be reproduced (or copied) under that name. Not significant -- it
    # identifies an input, not a product, and law's store paths are about products.
    hists_version = luigi.Parameter(
        default="",
        significant=False,
        description="version of the Hists_merged tree to read; defaults to --version",
    )

    @property
    def input_hists_version(self):
        return self.hists_version or self.version

    def output_dir_target(self, *path):
        """remote_dir_target() that also works when fs_default is a local directory.

        FLAF's remote_dir_target() handles a str fs and a remote fs but not a
        law.LocalFileSystem: it falls through to WLCGDirectoryTarget(), which raises
        "fs must be a RemoteFileSystem instance". Its sibling remote_target() (the file
        variant) already carries the branch mirrored here. The CI configures fs_default
        as a plain local path, so every directory output in this chain has to cope --
        without this, output() raises and luigi reports the task as having "error in
        complete() method".
        """
        fs = self.fs_default
        if isinstance(fs, law.LocalFileSystem):
            return law.LocalDirectoryTarget(os.path.join(*path), fs=fs)
        return self.remote_dir_target(*path)

    def datacard_config_path(self):
        return os.path.join(
            self.ana_path(), self.global_params["StatInference"]["config"]
        )

    def get_config_data(self):
        # Cached per instance: requires()/workflow_requires() are re-entered many times
        # during graph construction and each call would otherwise re-parse the yaml.
        if getattr(self, "_config_data", None) is None:
            with open(self.datacard_config_path(), "r") as f:
                self._config_data = yaml.safe_load(f)
        return self._config_data

    def get_era_groups(self):
        return self.get_config_data().get("era_groups", {})

    def get_top_level_eras(self):
        """Eras that carry their own datacards and limits: meta-eras plus any standalone
        real era. Real eras already covered by a meta-era are excluded so they are not
        double-counted alongside it."""
        data = self.get_config_data()
        eras = data.get("eras", [self.period])
        grouped = {e for sub in self.get_era_groups().values() for e in sub}
        return [e for e in eras if e not in grouped]

    def get_campaign(self, era):
        """dhi campaign key for an era, from the config's ``campaigns`` map.

        It cannot be derived from the era name. 'Run3_Early' happens to lowercase into
        dhi's 'run3_early', but 'Run3_2022' would give 'run3_2022', which is not a key of
        dhi.config.campaign_labels -- and dhi draws an unknown key verbatim rather than
        failing, so a guessed value ends up printed on the plot as the luminosity label.
        """
        return self.get_config_data().get("campaigns", {}).get(era)

    def get_masses(self):
        """Every parameter value the configuration declares, across all processes.

        Not restricted to signals: with ``param_dependent_bkg`` a background can be
        parameterised too, and its histograms live in the same per-mass input file.
        """
        masses = set()
        for proc in self.get_config_data().get("processes", []):
            if not isinstance(proc, dict):
                # A bare string entry names a background directly and carries no
                # settings -- the shape dc_make's load_config functions also accept.
                continue
            for value in proc.get("param_values") or []:
                masses.add(value)
        if not masses:
            raise RuntimeError(
                f"No process in {self.datacard_config_path()} declares param_values, so "
                "there is no mass point to build inputs for."
            )
        return sorted(masses)

    def get_required_variables(self):
        """The Hists_merged variables this chain reads: one per mass point (the 2D
        histogram HistRebinTask slices, or an already-1D shape for a configuration that
        does not rebin), derived from the datacard config's input_file_pattern.

        This is the *complete* input list -- there is no intersection with global.yaml's
        histTuple_flavor variable list anywhere downstream, so a mass point present here
        is read whether or not the active flavour happens to mention it.
        """
        data = self.get_config_data()
        model = data["model"]
        pattern = model["input_file_pattern"]
        param_name = model["parameters"][0]

        variables = set()
        for mass in self.get_masses():
            rel = Template(pattern).safe_substitute({"ERA": "ERA", param_name: mass})
            # "<era>/<variable>/<variable>.root": the variable is the directory holding
            # the file, not the file itself.
            variables.add(os.path.basename(os.path.dirname(rel)))
        return sorted(variables)

    def merged_hist_reqs(self, eras):
        """{(era, variable): MergedHists} -- every merged input file for those eras.

        Shared by the two tasks that read Hists_merged: HistRebinTask (over its discovery
        eras) and, for a configuration with no rebinning step, CreateDatacardsTask (over
        its sub-periods).
        """
        # Deferred: MergedHists subclasses this class, so importing it at module level
        # would be circular.
        from .MergedHists import MergedHists

        return {
            (era, variable): MergedHists.req(self, period=era, variable=variable)
            for era in eras
            for variable in self.get_required_variables()
        }

    @contextlib.contextmanager
    def stage_inputs(self, targets):
        """Assemble the required merged histograms into the "<era>/<variable>/
        <variable>.root" layout that Model.getInputFileName resolves input_file_pattern
        against, and yield the directory holding it.

        Each input is localized individually so only the files this task actually reads
        cross the network -- the merged tree holds every variable of the active
        histTuple_flavor (~106 for 'default') across every era.
        """
        with contextlib.ExitStack() as stack:
            staging = law.LocalDirectoryTarget(is_tmp=True)
            staging.touch()
            stack.callback(lambda: staging.remove(silent=True))

            for (era, variable), target in targets.items():
                dest_dir = os.path.join(staging.abspath, era, variable)
                os.makedirs(dest_dir, exist_ok=True)
                local_inp = stack.enter_context(target.localize("r"))
                shutil.copy2(
                    local_inp.abspath, os.path.join(dest_dir, f"{variable}.root")
                )

            yield staging

    def check_inputs(self, targets):
        """Report every missing input at once, naming the version they were looked for
        under. law already refuses to run with an incomplete MergedHists dependency; this
        makes the condor log self-explanatory when the task is forced anyway, instead of
        failing on whichever file ROOT happened to open first."""
        missing = sorted(
            f"{era}/{var}" for (era, var), t in targets.items() if not t.exists()
        )
        if missing:
            raise RuntimeError(
                f"{len(missing)} of {len(targets)} merged histograms are missing under "
                f"hists_version='{self.input_hists_version}' "
                f"(<fs_HistTuple>/{self.input_hists_version}/Hists_merged/<era>/<var>/<var>.root): "
                + ", ".join(missing)
            )
