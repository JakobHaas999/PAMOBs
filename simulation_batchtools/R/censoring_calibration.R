# Total-censoring calibration shared by null and level-shift runs.
# A population depends only on the data-generating parameters below. Sample size,
# spline df, repetition, seeds and fitting options do not change its calibration.
censoring_population_columns <- function() {
  c(
    "scenario", "beta0", "beta1", "beta2", "beta3", "shape", "rate",
    "meanlog", "sdlog", "covariate_setting", "gamma1", "gamm2",
    "delta", "tau", "admin_time"
  )
}

censoring_metadata <- function(design) {
  checkmate::assert_data_frame(design, min.rows = 1L)
  columns <- c(censoring_population_columns(), "censoring_rate")
  missing_columns <- setdiff(columns, names(design))
  if (length(missing_columns)) {
    stop("Missing calibration parameters: ", paste(missing_columns, collapse = ", "))
  }
  metadata <- as.data.frame(design)[columns]
  for (column in c("scenario", "covariate_setting")) {
    metadata[[column]] <- as.character(metadata[[column]])
    checkmate::assert_character(metadata[[column]], any.missing = FALSE)
  }
  for (column in setdiff(columns, c("scenario", "covariate_setting"))) {
    checkmate::assert_numeric(metadata[[column]], finite = TRUE)
  }
  checkmate::assert_numeric(metadata$delta, any.missing = FALSE, finite = TRUE)
  checkmate::assert_numeric(
    metadata$admin_time, any.missing = FALSE, finite = TRUE, lower = 1e-8
  )
  if (any(metadata$admin_time <= 1e-8)) {
    stop("admin_time must be greater than the numerical lower bound (1e-8).")
  }
  checkmate::assert_numeric(
    metadata$censoring_rate, any.missing = FALSE, lower = 0, upper = 1
  )
  if (any(metadata$censoring_rate >= 1)) {
    stop("Censoring targets must be smaller than one.")
  }
  shifted <- metadata$delta != 0
  if (any(shifted & (is.na(metadata$tau) | metadata$tau <= 0 |
                     metadata$tau >= metadata$admin_time))) {
    stop("A nonzero level shift requires 0 < tau < admin_time.")
  }
  metadata
}

# Integer codes preserve exact numeric values and give NA an explicit level.
# Compute keys jointly whenever rows from two different tables are compared.
censoring_row_keys <- function(metadata) {
  codes <- lapply(metadata, function(x) match(x, unique(x)))
  do.call(paste, c(codes, list(sep = ":")))
}

# make_input(pars, n) must return list(formula = ..., data = ...).
# The caller sources R/sim_surv.R before invoking this helper.
calibrate_study_censoring <- function(
  design, make_input, n_pilot = 50000L, seed = 900000L
) {
  checkmate::assert_function(make_input)
  checkmate::assert_count(n_pilot, positive = TRUE)
  checkmate::assert_count(seed)
  settings <- unique(censoring_metadata(design))
  population_keys <- censoring_row_keys(settings[censoring_population_columns()])
  groups <- split(
    seq_len(nrow(settings)),
    factor(population_keys, levels = unique(population_keys))
  )
  if (seed + length(groups) > .Machine$integer.max) {
    stop("Calibration seeds exceed the supported integer range.")
  }

  calibration <- lapply(seq_along(groups), function(i) {
    group <- settings[groups[[i]], , drop = FALSE]
    pars <- as.list(group[1L, , drop = FALSE])
    calibration_seed <- seed + i
    message(sprintf(
      "Calibrating population %d/%d: %s, %s, hazard factor %.3g",
      i, length(groups), pars$scenario, pars$covariate_setting, exp(pars$delta)
    ))
    set.seed(calibration_seed)
    input <- make_input(pars, n = n_pilot)
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
    if (any(group$censoring_rate < admin_fraction)) {
      stop(sprintf(
        paste0(
          "Censoring target %.4f is below administrative censoring %.4f ",
          "for %s, %s, delta %.4f, tau %.4f."
        ),
        min(group$censoring_rate), admin_fraction, pars$scenario,
        pars$covariate_setting, pars$delta, pars$tau
      ))
    }
    group$lambda_c <- vapply(group$censoring_rate, function(target) {
      calibrate_censoring(
        event_times = events, censoring_rate = target, admin_time = pars$admin_time
      )
    }, numeric(1))
    checkmate::assert_numeric(group$lambda_c, any.missing = FALSE, finite = TRUE, lower = 0)
    group$admin_fraction <- admin_fraction
    group$calibration_seed <- calibration_seed
    group
  })

  calibration <- do.call(rbind, calibration)
  rownames(calibration) <- NULL
  attr(calibration, "n_pilot") <- n_pilot
  attr(calibration, "seed") <- seed
  calibration
}

# Attach rates by exact population parameters and censoring target. In particular,
# changing df, sample size or setting_id does not invalidate an existing rate.
attach_study_calibration <- function(design, calibration) {
  metadata <- censoring_metadata(design)
  calibration_metadata <- censoring_metadata(calibration)
  if (!"lambda_c" %in% names(calibration)) {
    stop("Calibration metadata is incomplete: lambda_c is missing.")
  }
  checkmate::assert_numeric(
    calibration$lambda_c, any.missing = FALSE, finite = TRUE, lower = 0
  )
  keys <- censoring_row_keys(rbind(metadata, calibration_metadata))
  design_keys <- keys[seq_len(nrow(metadata))]
  calibration_keys <- keys[nrow(metadata) + seq_len(nrow(calibration_metadata))]

  first <- match(calibration_keys, calibration_keys)
  if (any(calibration$lambda_c != calibration$lambda_c[first])) {
    stop("Conflicting censoring rates for the same calibration parameters.")
  }
  index <- match(design_keys, calibration_keys)
  if (anyNA(index)) {
    stop(
      "Calibration does not match all design populations and censoring targets. ",
      "Use a new calibration file for the changed data-generating parameters."
    )
  }
  # Also explicitly check every metadata column without numerical tolerance.
  matches <- vapply(names(metadata), function(column) {
    isTRUE(all.equal(
      metadata[[column]], calibration_metadata[[column]][index],
      check.attributes = FALSE, tolerance = 0
    ))
  }, logical(1))
  if (!all(matches)) stop("Calibration metadata does not exactly match the design.")

  design$lambda_c <- calibration$lambda_c[index]
  design
}

# Loading a run can call this once: reuse a compatible cache or create it first.
# An existing incompatible cache is rejected and is never overwritten silently.
load_or_calibrate_censoring <- function(
  design, make_input, calibration_file, n_pilot = 50000L, seed = 900000L
) {
  checkmate::assert_string(calibration_file)
  checkmate::assert_count(n_pilot, positive = TRUE)
  checkmate::assert_count(seed)
  if (file.exists(calibration_file)) {
    calibration <- readRDS(calibration_file)
    if (!isTRUE(all.equal(attr(calibration, "n_pilot"), n_pilot, tolerance = 0))) {
      stop("The cached calibration uses a different or unknown n_pilot; use a new file.")
    }
    if (!isTRUE(all.equal(attr(calibration, "seed"), seed, tolerance = 0))) {
      stop("The cached calibration uses a different or unknown seed; use a new file.")
    }
    return(attach_study_calibration(design, calibration))
  }

  calibration <- calibrate_study_censoring(
    design = design, make_input = make_input, n_pilot = n_pilot, seed = seed
  )
  design <- attach_study_calibration(design, calibration)
  dir.create(dirname(calibration_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(calibration, calibration_file)
  design
}
