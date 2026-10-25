# Batchtools simulations

Run the scripts from the project root.

```text
simulation_batchtools/
  designs/
    pilot/pemtree_quadratic_pilot.R
    main/
  runs/
    pilot/pemtree_quadratic_pilot.R
    main/
  registries/
    pilot/
    main/
  results/
    pilot/
    main/
```

Each design defines `config$phase` (`"pilot"` or `"main"`) and a unique
`config$run_id`. Registries are created under `registries/<phase>/<run_id>/`.
Store collected results under `results/<phase>/<run_id>/`.

To register the current pilot study:

```r
source("simulation_batchtools/runs/pilot/pemtree_quadratic_pilot.R")
```

Run scripts create the registry when needed and submit the jobs. Registry data
and results are excluded from Git. The functions in `R/`, `dgps/`, `methods/`,
and `metrics/` are shared across both phases.

## Level-shift study

The shared design is `designs/level_shift_scenarios_base.R`. It combines the
three baseline hazards, two sample sizes, three censoring targets and two
covariate settings with hazard factors 1.25, 1.5, 2 and 3 (`delta = log(factor)`).
The change point is `tau = 5`, and administrative follow-up ends at 10.
This gives 144 settings. The pilot has 100 repetitions per setting (14,400 jobs);
main has 1,000 repetitions (144,000 jobs). Both use two workers by default,
configured in their respective design entry scripts.

First run the lightweight checks. They use temporary files and small simulations,
not the production registries or calibration:

```r
source("simulation_batchtools/tests/level_shift_smoke.R")
```

Prepare the censoring calibration once, then start the pilot:

```r
source("simulation_batchtools/calibration/level_shift_scenarios.R")
source("simulation_batchtools/runs/pilot/level_shift_scenarios_pilot.R")
```

Calibration uses 50,000 subjects per baseline/covariate/shift combination
(24 populations). Its RDS file is shared by pilot and main. It targets total
censoring, including administrative censoring; infeasible targets cause an
error. Calibration is never started implicitly when loading a design.
Metadata checks reject a calibration made for different design parameters.

After reviewing the pilot results, start main:

```r
source("simulation_batchtools/runs/main/level_shift_scenarios_main.R")
```

Keep the same `deltas` and `tau` in the shared design for calibration, pilot and
main. If these or other generating parameters change, regenerate calibration
and choose new run IDs for new studies; an existing registry retains its old
jobs. Set `config$run_id` after `make_level_shift_config()` in the relevant design
entry script.

The calibration RDS is excluded from Git. When moving the project to a server,
copy it separately or run the calibration there before starting a study.
