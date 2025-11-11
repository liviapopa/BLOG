library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# --- Files and target info ---
datasets <- list(
  "pos_median_mat30"  = list(file = "pos_median_mat_30_pi99.csv",  num_targets = 10),
  "pos_median_mat100" = list(file = "pos_median_mat_100_pi99.csv", num_targets = 20),
  "pos_median_mat350" = list(file = "pos_median_mat_350_pi99.csv", num_targets = 50)
)

# --- Function to compute metrics ---
compute_metrics <- function(column_data, num_targets) {
  sorted_data <- sort(abs(column_data), decreasing = TRUE)
  TP <- sum(sorted_data[1:num_targets] != 0)
  FP <- sum(sorted_data[(num_targets+1):length(column_data)] != 0)
  TN <- sum(sorted_data[(num_targets+1):length(column_data)] == 0)
  FN <- sum(sorted_data[1:num_targets] == 0)
  
  tibble(
    FPR = ifelse(FP + TN > 0, FP / (FP + TN), NA),
    TPR = ifelse(TP + FN > 0, TP / (TP + FN), NA),
    TNR = ifelse(TN + FP > 0, TN / (TN + FP), NA),
    FNR = ifelse(FN + TP > 0, FN / (FN + TP), NA),
    Precision = ifelse(TP + FP > 0, TP / (TP + FP), NA),
    Accuracy = (TP + TN) / (TP + TN + FP + FN)
  )
}

# --- Read files and compute metrics ---
all_metrics <- map2_dfr(
  names(datasets),
  datasets,
  ~{
    mat <- read.csv(.y$file)
    num_targets <- .y$num_targets
    
    map_dfr(seq_len(ncol(mat)), function(s) {
      compute_metrics(mat[, s], num_targets) %>% mutate(Dataset = .x, Simulation = s)
    })
  }
)

# --- Summarize metrics ---
summary_metrics <- all_metrics %>%
  group_by(Dataset) %>%
  summarise(across(FPR:Accuracy, list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE))))

print(summary_metrics)

# --- Generic violin plot function ---
plot_metric <- function(df, metric, ylims = NULL) {
  p <- ggplot(df, aes(x = Dataset, y = {{ metric }}, fill = Dataset)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.3) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white") +
    labs(x = "Simulation Scenario",
         y = deparse(substitute(metric)),
         title = paste0(deparse(substitute(metric)), " Across Simulations")) +
    theme_light() +
    theme(legend.position = "top") +
    scale_fill_manual(values = c("pos_median_mat30" = "blue",
                                 "pos_median_mat100" = "red",
                                 "pos_median_mat350" = "green"),
                      labels = c("30 Variables", "100 Variables", "350 Variables")) +
    scale_x_discrete(labels = c("pos_median_mat30" = "Posterior Medians 30 Variables",
                                "pos_median_mat100" = "Posterior Medians 100 Variables",
                                "pos_median_mat350" = "Posterior Medians 350 Variables"))
  
  if (!is.null(ylims)) p <- p + ylim(ylims)
  
  print(p)
}

# --- Plot all metrics ---
plot_metric(all_metrics, metric = FPR, ylims = c(0, 0.15))
plot_metric(all_metrics, metric = TPR)
plot_metric(all_metrics, metric = Precision)
plot_metric(all_metrics, metric = Accuracy)
