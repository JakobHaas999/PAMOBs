extract_scores_pem <- function(fit, ...) {
  checkmate::assert_class(fit, "glm")
  sandwich::estfun(fit, ...)
}

aggregate_scores_pem <- function(scores, ped, ...) {
  checkmate::assert_matrix(scores)
  checkmate::assert_data_frame(ped)
  checkmate::assert_true("id" %in% colnames(ped))

  rowsum(scores, group = ped[["id"]], ...)
}
