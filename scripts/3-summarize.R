# Nick Williams
# Research Biostatistician
# Department of Population Health Sciences
# Weill Cornell Medicine

# packages ----------------------------------------------------------------

library(data.table)
library(ggplot2)
library(patchwork)
library(purrr)

# global ------------------------------------------------------------------

source(here::here("R", "generation.R"))
source(here::here("R", "truth.R"))

res <- here::here("results")

alpha <- 0.05

plot_res <- function(data, scenario, y) {
  title <- if (y == "bias" || y == "rel_std_error") {
    glue::glue("Scenario {scenario}")
  } else {
    NULL
  }
  
  y_lab <- if (scenario == 1) {
    if (y == "bias") {
      "Bias"
    } else if (y == "n_bias") {
      expression(sqrt("n") %*% "Bias")
    } else if (y == "n_mse_bound") {
      expression("n" %*% "MSE" / "Eff. Bound")
    } else if (y == "rel_std_error") {
      "Relative Std. Error"
    } else {
      "Coverage"
    }
  }
  
  ggplot(data[type == scenario], aes_string(x = "n", y = y, group = "estimator")) + 
    geom_line(alpha = 0.8) + 
    geom_point(aes_string(shape = "estimator")) + 
    labs(y = y_lab, 
         title = title, 
         linetype = "Estimator")
}

lty <- c(0, 16, 2, 4)
names(lty) <- c("IPW", "substitution", "TMLE", "SDR")

theme_bias <- list(
  scale_shape_manual(breaks = c("IPW", "substitution", "TMLE", "SDR"),
                        labels = c("IPW", "Sub.", "TML", "SDR"), 
                        values = lty),
  theme_bw(base_size = 8),
  theme(plot.title = element_text(hjust = 0.5, size = 8), 
        plot.tag = element_text(size = 8))
)

theme_covr <- list(
  scale_shape_manual(breaks = c("TMLE", "SDR"),
                        labels = c("TML", "SDR"), 
                        values = lty),
  theme_bw(base_size = 8),
  theme(plot.title = element_text(hjust = 0.5, size = 8), 
        plot.tag = element_text(size = 8))
)

# data import -------------------------------------------------------------

fits <- readRDS(file.path(res, "fits-exported.rds"))
fits <- fits[std.error < 1 | is.na(std.error) & estimate < 1, ]

# results -----------------------------------------------------------------

# true value estimates
set.seed(5622)
true_vals <- true(1e5, sample(1e8, 1))
truth <- true_vals$true
bound <- true_vals$bound

results <- 
  fits[, .(mean = mean(estimate), 
           bias = abs(mean(estimate - truth)), 
           n_bias = abs(mean(sqrt(n) * (estimate - truth))), 
           rel_std_error = (mean(std.error) / sqrt(bound / n)), 
           coverage = mean(map2_lgl(conf.low, conf.high, ~ dplyr::between(truth, .x, .y))), 
           prop_coverage = mean(qnorm(alpha / 2) < (estimate - truth) / (sqrt(bound) / sqrt(n)) &
                                  (estimate - truth) / (sqrt(bound) / sqrt(n)) < qnorm(1 - alpha / 2)), 
           n_mse_bound = mean(n * (estimate - truth)^2 / bound)), 
       by = .(estimator, type, n)]

# figures -----------------------------------------------------------------

list(
  plot_res(results, 1, "bias"),
  plot_res(results, 2, "bias") + 
    coord_cartesian(ylim = c(0, 0.05)),
  plot_res(results, 3, "bias") + 
    coord_cartesian(ylim = c(0, 0.06)),
  plot_res(results, 4, "bias") + 
    coord_cartesian(ylim = c(0, 0.05)),
  plot_res(results, 1, "n_bias") + 
    coord_cartesian(ylim = c(0, 1.5)),
  plot_res(results, 2, "n_bias") + 
    coord_cartesian(ylim = c(0, 3.25)),
  plot_res(results, 3, "n_bias") + 
    coord_cartesian(ylim = c(0, 3.25)),
  plot_res(results, 4, "n_bias") + 
    coord_cartesian(ylim = c(0, 3.25)),
  plot_res(results, 1, "n_mse_bound"),
  plot_res(results, 2, "n_mse_bound"),
  plot_res(results, 3, "n_mse_bound"),
  plot_res(results, 4, "n_mse_bound")
) %>% 
  wrap_plots(ncol = 4, nrow = 3) + 
  plot_layout(guides = "collect") & 
  plot_annotation(tag_levels = 'A') & 
  theme_bias

ggsave("bias-plot.png", width = 9, height = 5, dpi = 600)

list(
  plot_res(results, 1, "rel_std_error") + 
    coord_cartesian(ylim = c(0.75, 1)),
  plot_res(results, 2, "rel_std_error") + 
    coord_cartesian(ylim = c(0.75, 1)),
  plot_res(results, 3, "rel_std_error") + 
    coord_cartesian(ylim = c(0.5, 1)),
  plot_res(results, 4, "rel_std_error")  + 
    coord_cartesian(ylim = c(0.55, 1)),
  plot_res(results, 1, "coverage") + 
    coord_cartesian(ylim = c(0.8, 1)),
  plot_res(results, 2, "coverage") + 
    coord_cartesian(ylim = c(0, 1)),
  plot_res(results, 3, "coverage") + 
    coord_cartesian(ylim = c(0, 1)),
  plot_res(results, 4, "coverage") + 
    coord_cartesian(ylim = c(0, 1))
) %>% 
  wrap_plots(ncol = 4, nrow = 2) + 
  plot_layout(guides = "collect") & 
  plot_annotation(tag_levels = 'A') & 
  theme_covr

ggsave("coverage-plot.png", width = 8.25, height = 3.25, dpi = 600)
