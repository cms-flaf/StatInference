# StatInference — instructions for Copilot code review

Datacard construction and statistical inference for the CMS analyses built on
[FLAF](https://github.com/cms-flaf/FLAF): it turns merged histograms into Combine datacards,
optimises binning, and runs limits and fits. Used as a submodule by HH_bbtautau and HH_bbWW.

**Read `FLAF/.github/copilot-instructions.md` first** for the shared rules on what a useful review
comment looks like here and what not to flag. The rule that documentation ships in the same PR
applies here too, and is restated below — it works differently for this repository, which has no
documentation site of its own.

## What costs the most here

This repository produces the **numbers that go in the paper**. A mistake does not crash: it
yields a limit that is plausible and wrong, and there is no downstream check that catches it. A
bug in FLAF wastes grid time; a bug here can be published.

So the questions worth asking of a diff are: *does this change a number?*, *does it change which
events or processes enter a bin?*, and *does it weaken a check that exists to catch a bad
histogram?*

## Invariants

### Safety checks exist because histograms go wrong quietly

- `resolveNegativeBins` (`common/tools.py`) rejects a shape with a **negative integral** before it
  ever looks at whether individual bins are negative-within-error. The per-process flags
  `allow_zero_integral`, `allow_negative_integral`, `allow_negative_bins_within_error` and
  `max_n_sigma_for_negative_bins` each switch off part of that protection.
- **Treat a new or widened `allow_*` flag as a finding** unless the PR explains why that
  process's shape is legitimately in that state. The honest fix is usually more statistics or
  dropping the process from the card, not disabling the guard. An amc@NLO sample with negative
  event weights and too few events is *not* a reason to allow a negative integral — it is a
  reason not to use it.
- A shape that silently becomes empty is worse than one that raises. Prefer failing loudly.

### Datacards

- `dc_make/maker.py` resolves processes from the card; a bare string entry means
  `Process(name, name)`, while a mapping carries `hist_name`, `is_signal`, `is_data`,
  `param_values` and the channel list. Adding an entry in the wrong form changes what is read
  from the input file, not just how it is labelled.
- The card's `eras:` list is its scope. Running a datacard task for an era the card does not
  declare fails on a missing histogram — that is the card's scope, not a bug to be patched around.
- A signal process with model parameters needs `param_values`; `param_dependent_bkg` changes how
  backgrounds are matched to a parameter point. Check that a change to one is reflected in the
  other.
- Histogram names are built from the card and must match what the analysis actually wrote. A
  rename on either side is a silent failure at datacard time.

### Binning optimisation

Changes under `bin_opt/` alter which events land in which bin, and therefore the sensitivity.
Look for: bin edges derived from the signal shape being applied to backgrounds, an optimisation
that can return fewer bins than requested without saying so, and any change to
`signalFractionForRelevantBins` semantics.

### Reproducibility

The output of this repository is compared across eras, channels and analyses. A change that makes
results depend on dictionary ordering, filesystem ordering or an unseeded random draw is a
finding even when the numbers look fine in one run.

## Documentation must ship with the change

A PR must update the documentation **in the same PR** whenever it changes anything a user can
observe: a task or its arguments, a card key, a CLI flag, an output location, a default, or the
meaning of an existing option.

**This repository has no documentation site of its own.** The statistical-inference chain is
documented in the analyses (`HH_bbtautau/docs/stat_inference.md`,
`HH_bbWW/docs/stat_inference.md`) and, where it is framework-wide, in `FLAF/docs/`. A
user-visible change here therefore needs a companion PR against the analysis or the framework;
flag the absence of one rather than accepting an undocumented change on the grounds that this
repository has nowhere to put it.

## Do not flag

- Missing unit tests for code that needs Combine, CMSSW or real histogram files.
- The `law/tasks.py` + script split, or shelling out to `combine`.
- Formatting — the `formatting-check` workflow settles it.
- Repetition between the per-analysis card configurations; they are kept explicit on purpose.

## Repository facts

Verified 2026-08-27; re-check before relying on any of it.

| | |
|---|---|
| Layout | `dc_make/` (datacard construction: `maker.py`, `create_datacards.py`, `process.py`, `model.py`, `uncertainty.py`, `binner.py`), `bin_opt/` (binning optimisation and limits), `common/` (`tools.py`, `param_parse.py`), `law/tasks.py`, `fit_val/`, `config/`, `legacy/` |
| Consumers | HH_bbtautau and HH_bbWW use this as a submodule; H_mumu does not |
| Cards | analysis-side, e.g. `HH_bbWW/config/Datacards/CI_card.yaml` |
| Tests | none in this repository |
| Workflows | `formatting-check`, `repo-sanity-checks`, `trigger-flaf-integration`. Formatting is checked automatically — do not comment on it |
| Integration test | Triggered by `@cms-flaf-bot please test`; the pipeline configuration lives in `cms-flaf/FLAF_ci` |
