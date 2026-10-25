# Batchtools simulations

Run all commands from the project root. Each file in `runs/` defines one complete
study: its settings, design, data generator, calibration and execution.
Shared functions live in `R/`; existing analyses remain in `analysis/`.

```text
simulation_batchtools/
  runs/
    pilot/
      null_scenarios_pilot.R
      level_shift_scenarios_pilot.R
      pemtree_quadratic_pilot.R
    main/
      null_scenario_main.R
      level_shift_scenarios_main.R
      pemtree_quadratic_main.R
  R/
  analysis/
  registries/<phase>/<run_id>/
  results/<phase>/<run_id>/
```

## Settings and paired df comparisons

Edit the settings near the top of the corresponding run function:
`dfs`, `n_rep`, `ncpus`, `run_id`, the scenario table and the factorial grid.
Use the same scientific settings and df grid in pilot and main.

Null and level-shift studies currently use `dfs <- c(4L, 6L, 8L)`.
These are unpenalized cubic B-spline basis dimensions in `splines::bs()`.
Each simulated sample is fitted with every df. `dgp_id` identifies the data
condition, `setting_id` additionally distinguishes df, and the seed depends on
`dgp_id` and repetition only. Seeds intentionally repeat across df. Pilot seeds
are a subset of main seeds for the same study configuration.

| Study | Data conditions | df values | Pilot repetitions | Pilot fits | Main repetitions | Main fits |
|---|---:|---:|---:|---:|---:|---:|
| Null | 36 | 3 | 5 | 540 | 1,000 | 108,000 |
| Level shift | 144 | 3 | 100 | 43,200 | 1,000 | 432,000 |

Both studies combine cubic, Weibull and lognormal baseline hazards with
`n = c(500L, 1000L)`, total-censoring targets 0.3/0.5/0.7 and conditions
without/with covariates. Follow-up ends at 10. The level-shift study additionally
crosses hazard factors 1.25/1.5/2/3 (`delta = log(factor)`), with a jump at 5.
Covariate effects are 0.3 for standard-normal x1 and 0.5 for Bernoulli(0.5) x2.

The older quadratic study retains its original settings and run IDs.
The new null and level-shift registries use names ending in `_df_v2`.
Choose a new run ID when changing an existing study. Resuming a registry with a
different saved job design is rejected; changes to function implementations also
require a new run ID.

## Inspect and start a run

Sourcing a run script only defines its entry function. It does not calibrate or
submit jobs. The `prepare` action builds the design without producing artifacts:

```r
source("simulation_batchtools/runs/pilot/null_scenarios_pilot.R")
study <- run_null_scenarios_pilot(action = "prepare")
head(study$design)
table(study$design$df)
```

From a terminal, preview either study:

```bash
Rscript simulation_batchtools/runs/pilot/null_scenarios_pilot.R prepare
Rscript simulation_batchtools/runs/pilot/level_shift_scenarios_pilot.R prepare
```

Start one pilot study:

```bash
Rscript simulation_batchtools/runs/pilot/null_scenarios_pilot.R
Rscript simulation_batchtools/runs/pilot/level_shift_scenarios_pilot.R
```

Equivalently, after sourcing its file, call `run_null_scenarios_pilot()` or
`run_level_shift_scenarios_pilot()`. These calls calibrate censoring, create the
registry, then run jobs with two workers by default. They block until completion.

After reviewing the pilot, start the corresponding main study:

```bash
Rscript simulation_batchtools/runs/main/null_scenario_main.R
Rscript simulation_batchtools/runs/main/level_shift_scenarios_main.R
```

## Censoring calibration

Calibration is computed by `calibrate_design()` whenever a run starts. The
function adds `lambda_c` to the design; no calibration files are loaded or saved.
To obtain a calibrated design without starting jobs, call the entry function
with `action = "calibrate"`, for example:

```r
source("simulation_batchtools/runs/pilot/null_scenarios_pilot.R")
study <- run_null_scenarios_pilot(action = "calibrate")
head(study$design)
```

The calibration uses 50,000 subjects for each population: six null populations
or 24 level-shift populations. It targets total censoring, including
administrative censoring. Rates are shared across all df values, sample sizes
and repetitions of the same population and censoring target.

Populations are grouped by `scenario`, `covariate_setting` and `delta`. This
assumes that baseline parameters, covariate effects, `tau` and `admin_time` are
fixed within each group, as in the current runs. If these parameters are varied,
include them in the grouping as well.

Pilot and main recompute calibration using the same population order, seed and
calibration sample size. With the same generating settings and R environment,
they therefore obtain the same rates. No separate calibration file needs to be
transferred when moving a run to the server.

## Progress on a server

Use a persistent terminal session, for example `tmux`, for long Socket runs.
From a second R session in the same project directory:

```r
reg_view <- batchtools::loadRegistry(
  file.dir = paste0(
    "simulation_batchtools/registries/main/",
    "level_shift_scenario_main_df_v2"
  ),
  writeable = FALSE, conf.file = NA, make.default = FALSE
)
batchtools::getStatus(reg = reg_view)
done <- batchtools::findDone(reg = reg_view)
if (nrow(done)) {
  print(batchtools::reduceResultsList(ids = tail(done, 10L), reg = reg_view))
}
```

Reload the registry to refresh the status. The registry exists only after
calibration and job registration. Model errors caught by `run_one()` are stored
as `success = FALSE` and `error` in results; batchtools can still mark those jobs
as done. The runner submits jobs that have not yet been submitted. Interrupted
jobs already marked as submitted are not automatically resubmitted.

Registries and collected results are excluded from Git.
Transfer them separately when moving results between the server and your computer.

Note: AI generated