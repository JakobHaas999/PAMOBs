# Calibrate total censoring for each population in the level-shift design.
# A single event-time sample is shared across sample sizes and censoring targets.
# The returned table retains setting metadata so stale calibrations can be rejected.
level_shift_settings <- function(design) {
  checkmate::assert_data_frame(design, min.rows = 1L)
  if (!"setting_id" %in% names(design) || anyNA(design$setting_id)) {
    stop("The design must contain nonmissing setting_id values.")
  }

  metadata <- setdiff(names(design), c("rep", "seed", "lambda_c"))
  settings <- unique(design[metadata])
  if (anyDuplicated(settings$setting_id)) {
    stop("Each setting_id must identify exactly one parameter combination.")
  }
  settings <- settings[order(settings$setting_id), , drop = FALSE]
  rownames(settings) <- NULL
  settings
}

calibrate_level_shift <- function(design, n_pilot = 50000L, seed = 1900000L) {
  checkmate::assert_count(n_pilot, positive = TRUE)
  checkmate::assert_count(seed)
  settings <- level_shift_settings(design)

  population_columns <- c(
    "scenario_id", "scenario", "beta0", "beta1", "beta2", "beta3",
    "shape", "rate", "meanlog", "sdlog", "covariate_setting",
    "gamma1", "gamm2", "delta", "tau", "admin_time"
  )
  required <- c(population_columns, "censoring_rate")
  missing_columns <- setdiff(required, names(settings))
  if (length(missing_columns)) {
    stop("Missing calibration parameters: ", paste(missing_columns, collapse = ", "))
  }
  checkmate::assert_numeric(settings$censoring_rate, any.missing = FALSE, lower = 0)
  if (any(settings$censoring_rate >= 1)) {
    stop("Censoring targets must be smaller than one.")
  }
  checkmate::assert_numeric(settings$admin_time, any.missing = FALSE, finite = TRUE)
  if (any(settings$admin_time <= 1e-8)) {
    stop("admin_time must be greater than the numerical lower bound (1e-8).")
  }

  # match() gives unused baseline parameters (NA) an explicit group level.
  group_levels <- lapply(settings[population_columns], function(x) match(x, unique(x)))
  group_id <- do.call(interaction, c(group_levels, list(drop = TRUE, lex.order = TRUE)))
  groups <- split(seq_len(nrow(settings)), group_id)
  if (seed + length(groups) > .Machine$integer.max) {
    stop("Calibration seeds exceed the supported integer range.")
  }

  calibration <- lapply(seq_along(groups), function(i) {
    group <- settings[groups[[i]], , drop = FALSE]
    pars <- as.list(group[1L, , drop = FALSE])
    calibration_seed <- seed + i
    message(sprintf(
      "Calibrating population %d/%d: %s, %s, hazard factor %.3g, tau %.3g",
      i, length(groups), pars$scenario, pars$covariate_setting,
      exp(pars$delta), pars$tau
    ))
    set.seed(calibration_seed)
    input <- make_level_shift_input(pars, n = n_pilot)
    events <- sim_event_times(
      formula = input$formula,
      data = input$data,
      breaks = if (pars$delta == 0) numeric(0) else pars$tau,
      maxt = pars$admin_time,
      interval = c(1e-8, 1.01 * pars$admin_time)
    )
    if (
      !is.data.frame(events) ||
      !all(c("eventtime", "status") %in% names(events)) ||
      nrow(events) != n_pilot ||
      anyNA(events$eventtime) || any(!is.finite(events$eventtime)) ||
      any(events$eventtime < 0) || anyNA(events$status) ||
      any(!events$status %in% c(0L, 1L))
    ) {
      stop("Invalid event-time sample during censoring calibration.")
    }

    admin_fraction <- mean(events$eventtime >= pars$admin_time)
    targets <- unique(group$censoring_rate)
    if (any(targets < admin_fraction)) {
      stop(sprintf(
        paste0(
          "Censoring target %.4f is below administrative censoring %.4f ",
          "for %s, %s, delta %.4f, tau %.4f."
        ),
        min(targets), admin_fraction, pars$scenario,
        pars$covariate_setting, pars$delta, pars$tau
      ))
    }
    rates <- vapply(targets, function(target) {
      calibrate_censoring(
        event_times = events,
        censoring_rate = target,
        admin_time = pars$admin_time
      )
    }, numeric(1))
    if (anyNA(rates) || any(!is.finite(rates)) || any(rates < 0)) {
      stop("Calibration produced invalid censoring rates.")
    }

    group$lambda_c <- rates[match(group$censoring_rate, targets)]
    group$admin_fraction <- admin_fraction
    group$calibration_seed <- calibration_seed
    group
  })

  calibration <- do.call(rbind, calibration)
  calibration <- calibration[order(calibration$setting_id), , drop = FALSE]
  rownames(calibration) <- NULL
  attr(calibration, "n_pilot") <- n_pilot
  calibration
}

# Attach an existing calibration only when every setting's metadata agrees.
# This function never starts a calibration as a side effect of loading a design.
attach_level_shift_calibration <- function(design, calibration_file) {
  checkmate::assert_string(calibration_file)
  if (!file.exists(calibration_file)) {
    stop(
      "Calibration file does not exist: ", calibration_file,
      ". Run simulation_batchtools/calibration/level_shift_scenarios.R first."
    )
  }
  settings <- level_shift_settings(design)
  calibration <- readRDS(calibration_file)
  if (!is.data.frame(calibration)) {
    stop("The calibration file must contain a data frame.")
  }
  required <- c(names(settings), "lambda_c")
  missing_columns <- setdiff(required, names(calibration))
  if (length(missing_columns)) {
    stop(
      "Calibration metadata is incomplete: ", paste(missing_columns, collapse = ", "),
      ". Recompute the level-shift calibration."
    )
  }
  if (anyNA(calibration$setting_id) || anyDuplicated(calibration$setting_id)) {
    stop("Calibration setting_id values must be nonmissing and unique.")
  }
  checkmate::assert_numeric(
    calibration$lambda_c, any.missing = FALSE, finite = TRUE, lower = 0
  )
  index <- match(settings$setting_id, calibration$setting_id)
  if (anyNA(index)) {
    stop("Calibration is missing design settings. Recompute the level-shift calibration.")
  }
  mismatched <- names(settings)[!vapply(names(settings), function(column) {
    isTRUE(all.equal(
      settings[[column]], calibration[[column]][index],
      check.attributes = FALSE, tolerance = 0
    ))
  }, logical(1))]
  if (length(mismatched)) {
    stop(
      "Calibration does not match the design (", paste(mismatched, collapse = ", "),
      "). Recompute the level-shift calibration."
    )
  }

  design$lambda_c <- calibration$lambda_c[match(design$setting_id, calibration$setting_id)]
  design
}
