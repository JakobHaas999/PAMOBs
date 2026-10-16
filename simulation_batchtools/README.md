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

This registers the jobs; it does not start them. Registry data and results are
excluded from Git. The functions in `R/`, `dgps/`, `methods/`, and `metrics/`
are shared across both phases.
