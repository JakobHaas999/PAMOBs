calibrate_design <- function(
  design,
  make_input,
  n_calibration = 50000L,
  seed = 900000L
) {
  # Gleiche Population: gleiche Baseline, Kovariaten und Sprunggröße.
  population <- interaction(
    design$scenario,
    design$covariate_setting,
    design$delta,
    drop = TRUE
  )

  groups <- split(
    seq_len(nrow(design)),
    factor(population, levels = unique(population))
  )

  design$lambda_c <- NA_real_

  for (i in seq_along(groups)) {
    rows <- groups[[i]]
    pars <- as.list(design[rows[1L], , drop = FALSE])

    message("Calibrating population ", i, "/", length(groups))
    set.seed(seed + i)

    # Große Stichprobe für diese Population simulieren.
    input <- make_input(pars, n = n_calibration)

    events <- sim_event_times(
      formula = input$formula,
      data = input$data,
      breaks = if (pars$delta == 0) numeric(0) else pars$tau,
      maxt = pars$admin_time,
      interval = c(1e-8, 1.01 * pars$admin_time)
    )

    # Für jedes Zensierungsziel eine Rate bestimmen.
    targets <- unique(design$censoring_rate[rows])

    rates <- vapply(
      targets,
      function(target) {
        calibrate_censoring(
          event_times = events,
          censoring_rate = target,
          admin_time = pars$admin_time
        )
      },
      numeric(1)
    )

    # Raten allen passenden Jobs zuordnen.
    design$lambda_c[rows] <- rates[
      match(design$censoring_rate[rows], targets)
    ]
  }

  design
}
