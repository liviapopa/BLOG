#========================
# Generalized Zellner g-function pipeline
#========================

library(MASS)

# Path to dataset (CSV or RData)
data_path <- "yourdataset.RData"  # or "yourdataset.csv"

# Name of the response variable in your dataset
response_var <- "Y"  

# Output file for Bayes factor ranking
output_file <- "bayes_factor_ranking.csv"

#------------------------
# Load dataset
#------------------------
if (grepl("\\.RData$", data_path)) {
  load(data_path)
  # Expect a data frame named 'data' or adjust accordingly
} else if (grepl("\\.csv$", data_path)) {
  data <- read.csv(data_path, header = TRUE, stringsAsFactors = FALSE)
} else {
  stop("Unsupported file type. Use .RData or .csv")
}

#------------------------
# Zellner g-prior function
#------------------------

zell_g <- function(Y, X, b0, g = NULL) {
  n <- length(Y)
  
  
  beta_hat <- solve(t(X) %*% X) %*% (t(X) %*% Y)
  Y_hat <- X %*% beta_hat
  Y0 <- X %*% b0
  sigmasq_hat <- norm(Y - X %*% beta_hat, type = "2")^2 / (n - p)
  
  if (is.null(g)) {
    g <- sqrt(n)
  } else if (g == "SURE") {
    g <- (norm(Y_hat - Y0, type = "2")^2 / (p * sigmasq_hat)) - 1
  }
  
  V0 <- g * solve(t(X) %*% X)
  b_star <- (1 / (1 + g)) * b0 + (g / (1 + g)) * beta_hat
  V_star <- (g / (1 + g)) * solve(t(X) %*% X)
  md <- t(b_star) %*% solve(V_star) %*% b_star
  linear_model <- lm(Y~0+X)
  r_sq <- summary(linear_model)$r.squared
  bf <- (1+g)^((n-p-1)/2)*(1+g*(1-r_sq))^(-(n-1)/2)
  
  return(list(mean = b_star, sd = V_star * sigmasq_hat, g = g,bf=bf))
}

#------------------------
# Prepare data for analysis
#------------------------
# Here we assume 'data' is a matrix-like object
# Users can modify this part for their specific dataset
dat_mat <- as.matrix(data)  # Convert to numeric matrix if needed
n_metabolites <- nrow(dat_mat)
metab_names <- rownames(dat_mat)

ranking_df <- data.frame(
  g_sqrtn = numeric(0),
  bayes_factor_sqrtn = numeric(0),
  g_sure = numeric(0),
  bayes_factor_sure = numeric(0),
  linear_r_sq_original_data = numeric(0)
)

for (metab in 1:dim(dat_array)[1]) {
  n <- ncol(dat_array[metab,,])   # number of samples
  tt <- nrow(dat_array[metab,,])  # number of timepoints
  
  # Y construction
  y_tilde <- t(dat_array[metab,,])
  DY_tilde <- array(0, c(n, tt - 1))
  for (i in 1:(tt - 1)) {
    DY_tilde[,i] <- y_tilde[,i + 1] - y_tilde[,i]
  }
  Y <- as.vector(DY_tilde)
  
  # X construction
  x_tilde <- dat_array[metab,,]
  DX_tilde <- array(0, c(n, tt - 1))
  for (i in 1:(tt - 1)) {
    DX_tilde[,i] <- x_tilde[,i + 1] - x_tilde[,i]
  }
  
  # Build block X matrix
  X <- array(0, c(n * (tt - 1), tt * (tt - 1) / 2))
  for (i in 1:ncol(DX_tilde)) {
    row_spot <- ((i - 1) * n) + 1
    col_spot <- 1 + (i * (i - 1) / 2)
    group <- matrix(nrow = n, ncol = 1)
    for (j in 1:i) {
      group <- cbind(group, DX_tilde[,j])
    }
    group <- group[,-1]
    X[row_spot:(row_spot + n - 1), col_spot:(col_spot + i - 1)] <- group
  }
  
  p <- ncol(X)
  
  # Zellner g-prior calculations
  zell_g_sqrtn <- zell_g(Y, X, b0 = rep(0, p))
  zell_g_sure <- zell_g(Y, X, b0 = rep(0, p), g = "SURE")
  
  # Original OLS R^2
  linear_model_original <- lm(Y ~ 0 + X)
  linear_r_sq_original <- summary(linear_model_original)$r.squared
  
  # Save results
  new_row <- data.frame(
    g_sqrtn = round(zell_g_sqrtn$g, 3),
    bayes_factor_sqrtn = round(zell_g_sqrtn$bf, 3),
    g_sure = round(zell_g_sure$g, 3),
    bayes_factor_sure = round(zell_g_sure$bf, 3),
    linear_r_sq_original_data = linear_r_sq_original
  )
  ranking_df <- rbind(ranking_df, new_row)
}


rownames(ranking_df) <- metab_names
View(ranking_df)

write.csv(ranking_df,"metabolites_bayes_factors.csv")
