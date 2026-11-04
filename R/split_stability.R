# Subsampling Stability of the Root Split
#
# Assesses the uncertainty of the root split of a model-based tree
# (e.g. pemtree, pamtree). The tree is fitted once to the full data and
# then refitted to B subsamples of individuals drawn without
# replacement. For each subsample, the function records whether the
# root node is split, on which variable, at which split point, and the
# p-value of the root parameter instability test.
#
# Subsampling is done on the level of individuals (rows of the original
# data), not on the PED rows, so that all intervals of an individual
# stay together. fit_fun rebuilds the PED for every subsample. Sampling
# without replacement avoids duplicated individuals, which the
# instability tests would treat as independent observations, making
# splits too likely. Fixed cut points must be supplied in fit_args$cut
# so that split points are comparable across subsamples.
#
# maxdepth is set to 2, so every tree has at most one split, and only
# the root split is evaluated.
#
# Input:
#   fit_fun   Tree fitting function, e.g. pemtree or pamtree. Must accept
#             a data argument and return a list with element tree
#             (a partykit modelparty).
#   fit_args  Named list of arguments for fit_fun. Must contain data
#             (one row per individual) and cut (fixed interval
#             boundaries). A supplied maxdepth is overwritten.
#   B         Number of subsamples (default: 500).
#   frac      Fraction of individuals per subsample, between 0.1 and
#             0.99 (default: 0.632).
#
# Output:
#   A list with
#     full_fit     Tree fitted to the full data.
#     split        Logical vector of length B: was the root node split?
#                  NA if the fit failed or the root model could not be
#                  estimated (e.g. too few events).
#     split_var    Name of the split variable; NA if there is no split
#                  or split is NA.
#     split_point  Split point for numeric split variables; NA if there
#                  is no split, split is NA, or the split is on a factor.
#     p_value      p-value of the root test as reported by
#                  partykit::info_node(); NA if not available.
#     idx          List of length B with the row indices of
#                  fit_args$data drawn in each subsample.
split_stability <- function(fit_fun, fit_args, B = 500, frac = 0.632) {
  checkmate::assert_function(fit_fun)
  checkmate::assert_list(fit_args)
  checkmate::assert_data_frame(fit_args$data)
  checkmate::assert_count(B, positive = TRUE)
  checkmate::assert_number(frac, lower = 0.1, upper = 0.99)
  if (is.null(fit_args$cut)) stop("Supply fixed cut points in fit_args$cut.")

  fit_args$maxdepth <- 2
  full_fit <- do.call(fit_fun, fit_args)

  D <- fit_args$data
  fit_args$data <- NULL
  N <- nrow(D)
  subsample_size <- floor(frac * N)

  split <- logical(B)
  split_point <- numeric(B)
  split_var <- character(B)
  p_value <- rep(NA_real_, B)
  idx_list <- vector("list", length = B)

  for (b in seq_len(B)) {
    # Subsample data
    idx <- sample.int(N, size = subsample_size, replace = FALSE)
    idx_list[[b]] <- idx
    D_star <- D[idx, , drop = FALSE]

    # Fit tree with subsampled data
    fit_star <- tryCatch(do.call(fit_fun, c(list(data = D_star), fit_args)), error = function(e) NULL)

    if (is.null(fit_star)) {
      split[b] <- NA
      split_point[b] <- NA_real_
      split_var[b] <- NA_character_
      next
    }

    tree_star <- fit_star$tree
    root <- partykit::node_party(tree_star)
    info <- partykit::info_node(root)
    if (!is.finite(info$objfun)) {
      split[b] <- NA
      split_point[b] <- NA_real_
      split_var[b] <- NA_character_
      next
    }
    if (length(info$p.value) == 1) p_value[b] <- info$p.value

    split[b] <- !partykit::is.terminal(root)
    if (split[b]) {
      sn <- partykit::split_node(root)
      split_var[b] <- names(tree_star$data)[partykit::varid_split(sn)]
      brk <- partykit::breaks_split(sn)
      split_point[b] <- if (is.null(brk)) NA_real_ else brk
    } else {
      split_point[b] <- NA_real_
      split_var[b] <- NA_character_
    }
  }

  list(
    full_fit = full_fit,
    split = split,
    split_var = split_var,
    split_point = split_point,
    p_value = p_value,
    idx = idx_list
  )
}
