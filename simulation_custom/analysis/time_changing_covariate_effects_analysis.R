source("setup.R")
## Helper functions --------------------------------------
get_splits <- function(tree) {
  inner_ids <- setdiff(nodeids(tree), nodeids(tree, terminal = TRUE))
  if (!length(inner_ids)) {
    return(numeric(0))
  }
  sort(unlist(nodeapply(
    tree,
    ids = inner_ids,
    FUN = function(node) breaks_split(split_node(node))
  ), use.names = FALSE))
}

## Analysis of results -----------------------------------
results_dir <- "simulation_custom/results"

## Scenario 1
results_1 <- readRDS(file.path(results_dir, "time_changing_effects_1.rds"))
trees_1 <- results_1$trees
hazards_1 <- results_1$hazards

splits_1 <- lapply(trees_1, get_splits)
mean(length(splits_1) > 1)
summary(unlist(splits_1))

# Plot hazards
hazards_1 <- as.data.table(hazards_1)
hazards_1 <- hazards_1[, .(
  mean = mean(hazard),
  lower = quantile(hazard, 0.025),
  upper = quantile(hazard, 0.975)
), by = c("tend", "x1", "x2")]
hazards_1[, let(
  truth = exp(log_hazard0(tend) + 0.4 * x1 + 0.6 * x2 + 0.8 * x1 * (tend >= 4)),
  x1 = factor(x1)
)]

Set1 <- RColorBrewer::brewer.pal(9, "Set1")
ggplot(hazards_1, aes(x = tend, col = x1, fill = x1)) +
  geom_stepribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, col = NA) +
  geom_stephazard(aes(y = mean)) +
  geom_line(aes(y = truth), lty = 2) +
  geom_vline(xintercept = tau, lty = 3) +
  scale_colour_manual(values = Set1[c(7, 5, 4)]) +
  scale_fill_manual(values = Set1[c(7, 5, 4)]) +
  labs(
    x = expression(t), col = expression(x[1]), fill = expression(x[1]),
    y = expression(lambda(t))
  ) +
  theme_bw()
