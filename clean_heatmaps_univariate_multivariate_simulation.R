library(reshape2)
library(ggplot2)
library(grid)

# ----------------------------
# Load all CSVs into a list
# ----------------------------
file_list <- list(
  sure_dat30 = "zell_bf_sure_mat_30.csv",
  sure_dat100 = "zell_bf_sure_mat_100.csv",
  sure_dat350 = "zell_bf_sure_mat_350.csv",
  sqrt_dat30 = "zell_bf_sqrt_mat_30.csv",
  sqrt_dat100 = "zell_bf_sqrt_mat_100.csv",
  sqrt_dat350 = "zell_bf_sqrt_mat_350.csv",
  gll_count_30 = "gllcount_30_pi99.csv",
  gll_count_100 = "gllcount_100_pi99.csv",
  gll_count_350 = "gllcount_350_pi99.csv",
  pos_median30 = "pos_median30_pi99.csv",
  pos_median100 = "pos_median100_pi99.csv",
  pos_median350 = "pos_median350_pi99.csv"
)

data <- lapply(file_list, read.csv)

# Transpose matrices where needed
data$gll_count_350 <- t(data$gll_count_350)
data$pos_median350 <- t(data$pos_median350)

# ----------------------------
# Function to clean a matrix
# ----------------------------
clean_matrix <- function(mat) {
  mat <- as.matrix(mat)
  mat[is.na(mat) | is.nan(mat) | is.infinite(mat)] <- 0
  return(mat)
}

# ----------------------------
# Generalized Heatmap Function
# ----------------------------
plot_heatmap <- function(mat, title = "Heatmap", x_breaks = NULL) {
  mat <- clean_matrix(mat)
  df_melt <- melt(mat)
  
  df_melt$color_category <- cut(
    df_melt$value,
    breaks = c(-Inf, 99, 300, 2000, 14902, Inf),
    labels = c("0 to 99", "99 to 300", "300 to 2000", "2000 to 14902", "above 14902")
  )
  
  color_values <- c(
    "0 to 99" = "white",
    "99 to 300" = "lightpink",
    "300 to 2000" = "hotpink",
    "2000 to 14902" = "red",
    "above 14902" = "darkred"
  )
  
  p <- ggplot(df_melt, aes(x = Var2, y = Var1, fill = color_category)) +
    geom_tile() +
    scale_fill_manual(
      values = color_values,
      name = "Value Range",
      labels = c(
        "white: 0 to 99", "pink: 99 to 300",
        "hot pink: 300 to 2000", "red: 2000 to 14902",
        "dark red: above 14902"
      )
    ) +
    labs(x = "Simulations", y = "Number of Variables", title = title) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
  
  if (!is.null(x_breaks)) {
    p <- p + scale_x_discrete(breaks = x_breaks)
  }
  
  print(p)
}

# ----------------------------
# Change only this line to choose which heatmap to plot
# ----------------------------
selected_dataset <- "sure_dat350"

# Display heatmap
plot_heatmap(data[[selected_dataset]], title = paste0("Heatmap", selected_dataset), x_breaks = c("1","25","50","100"))

# ----------------------------
# Generic function to plot single-row heatmap
# ----------------------------
plot_heatmap_row <- function(values, main_title) {
  values <- as.numeric(unlist(values))
  matrix_data <- matrix(values, nrow = 1)
  normalized_values <- (values - min(values)) / (max(values) - min(values))
  color_palette <- colorRampPalette(c("white", "darkred"))(100)
  color_indices <- round(normalized_values * 99) + 1
  colors <- color_palette[color_indices]
  
  plot(1:ncol(matrix_data), rep(1, ncol(matrix_data)), type = "n",
       xlab = "", ylab = "", axes = FALSE, ylim = c(0, 1.2))
  
  rect_height <- 0.8
  for (i in 1:ncol(matrix_data)) {
    rect(i - 0.5, 0, i + 0.5, rect_height, col = colors[i], border = "black")
  }
  
  title(main = main_title, line = 0.5)
  axis(1, at = 1:ncol(matrix_data), labels = 1:ncol(matrix_data), cex.axis = 0.7, las = 1)
  mtext("Number of Variables", side = 1, line = 2.5)
}

plot_single_heatmap <- function(selected_values, main_title) {
  plot_heatmap_row(selected_values, main_title)
}

# ----------------------------
# Change this to choose which dataset to plot
# ----------------------------
selected_values <- data$gll_count_350
main_title <- "Coefficients (350 Variables)"

plot_single_heatmap(selected_values, main_title)

# ----------------------------
# Plot multiple datasets
# ----------------------------
plot_multiple_heatmaps <- function(list_of_values, titles) {
  par(mfrow = c(length(list_of_values), 1), mar = c(4, 4, 2, 1))
  for (i in seq_along(list_of_values)) {
    plot_heatmap_row(list_of_values[[i]], titles[i])
  }
}

# Posterior medians
plot_multiple_heatmaps(
  list(data$pos_median30, data$pos_median100, data$pos_median350),
  c("Posterior Medians (30 Variables)", 
    "Posterior Medians (100 Variables)", 
    "Posterior Medians (350 Variables)")
)

# Coefficients
plot_multiple_heatmaps(
  list(data$gll_count_30, data$gll_count_100, data$gll_count_350),
  c("Coefficients (30 Variables)", 
    "Coefficients (100 Variables)", 
    "Coefficients (350 Variables)")
)
