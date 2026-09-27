sim_pexp_own <- function(formula, data, cut) {
  data <- data |> dplyr::mutate(
    id = dplyr::row_number(), time = max(cut),
    status = 1
  )

  Form <- Formula::Formula(formula)
  f1 <- stats::formula(Form, rhs = 1)

  if (length(Form)[2] > 1L) {
    stop("Two-part formulas are not supported by sim_pexp_own().")
  }

  ped <- pammtools::split_data(
    formula = survival::Surv(time, status) ~ .,
    data = dplyr::select_if(data, is.atomic),
    cut = cut,
    id = "id"
  ) |>
    dplyr::rename(t = "tstart")

  eta <- lazyeval::f_eval(f1, ped)
  pammtools:::assert_numeric_eta(eta, f1, data)
  ped[["rate"]] <- exp(eta)

  sim_df <- ped |>
    dplyr::group_by(id) |>
    dplyr::summarise(
      time = pammtools:::rpexp(
        rate = .data$rate,
        t = .data$t
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(status = 1L)

  sim_df <- suppressMessages(
    dplyr::left_join(
      sim_df,
      dplyr::select(
        data,
        -dplyr::all_of(c("time", "status"))
      ),
      by = "id"
    )
  )

  as.data.frame(sim_df)
}
