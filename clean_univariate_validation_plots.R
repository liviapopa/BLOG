library(dplyr)
library(tibble)
library(ggplot2)
library(cowplot)
library(purrr)

# --- Settings ---
threshold <- 99

# --- Files and targets ---
files <- list(
  "sqrt30" = "zell_bf_sqrt_mat_30.csv",
  "sure30" = "zell_bf_sure_mat_30.csv",
  "sqrt100" = "zell_bf_sqrt_mat_100.csv",
  "sure100" = "zell_bf_sure_mat_100.csv",
  "sqrt350" = "zell_bf_sqrt_mat_350.csv",
  "sure350" = "zell_bf_sure_mat_350.csv"
)

targets_per_scenario <- c(
  "sqrt30" = 10,
  "sure30" = 10,
  "sqrt100" = 20,
  "sure100" = 20,
  "sqrt350" = 50,
  "sure350" = 50
)

# --- Helper: Compute metrics ---
compute_metrics <- function(pred, truth) {
  TP <- sum(pred == 1 & truth == 1)
  FP <- sum(pred == 1 & truth == 0)
  TN <- sum(pred == 0 & truth == 0)
  FN <- sum(pred == 0 & truth == 1)
  
  FPR <- if ((FP + TN) == 0) NA else FP / (FP + TN)
  TPR <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  Precision <- if ((TP + FP) == 0) NA else TP / (TP + FP)
  Accuracy <- (TP + TN) / (TP + TN + FP + FN)
  
  tibble(TP, FP, TN, FN, FPR, TPR, Precision, Accuracy)
}

# --- Load data and compute metrics ---
all_results <- map2(names(files), files, ~{
  dat <- read.csv(.y, header = TRUE, check.names = FALSE)
  n_targets <- targets_per_scenario[.x]
  truth <- c(rep(1, n_targets), rep(0, nrow(dat) - n_targets))
  
  map_dfr(seq_len(ncol(dat)), function(s) {
    bf_col <- as.numeric(dat[, s])
    keep <- !is.na(bf_col)
    pred <- ifelse(bf_col[keep] > threshold, 1, 0)
    
    tibble(
      scenario = .x,
      simulation = s
    ) %>%
      bind_cols(compute_metrics(pred, truth[keep]))
  })
})
results_df <- bind_rows(all_results) %>%
  mutate(
    scenario = factor(scenario, levels = names(files)),
    type = ifelse(grepl("^sqrt", scenario), "sqrt n", "sure")
  )

# --- Summaries per scenario ---
summaries <- results_df %>%
  group_by(scenario) %>%
  summarise(
    n_sims = n(),
    across(c(FPR, TPR, Precision, Accuracy),
           list(
             mean = ~mean(.x, na.rm = TRUE),
             sd = ~sd(.x, na.rm = TRUE),
             median = ~median(.x, na.rm = TRUE),
             IQR_low = ~quantile(.x, 0.25, na.rm = TRUE),
             IQR_high = ~quantile(.x, 0.75, na.rm = TRUE)
           ),
           .names = "{.col}_{.fn}")
  )

print(summaries)

# --- Helper: Plot a metric ---
plot_metric <- function(df, y_var, y_label, ylim = NULL) {
  ggplot(df, aes(x = scenario, y = .data[[y_var]], fill = type)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    stat_summary(fun = "mean", geom = "point", shape = 23, size = 3, fill = "white") +
    scale_fill_manual(name = "Choice of g", values = c("sqrt n" = "blue", "sure" = "red")) +
    labs(title = paste(y_label, "(BF >", threshold, ")"), y = y_label, x = "Simulation Scenario") +
    ylim(ylim) +
    theme_minimal()
}

# --- Generate plots ---
p_fpr  <- plot_metric(results_df, "FPR", "False Positive Rate", c(0, NA))
p_prec <- plot_metric(results_df, "Precision", "Precision", c(0.6, 1))
p_acc  <- plot_metric(results_df, "Accuracy", "Accuracy", c(0.6, 1))

plot_grid(p_fpr, p_prec, p_acc, ncol = 1)
