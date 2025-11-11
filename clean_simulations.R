# ----------------------------
# Libraries
# ----------------------------
library(MBSGS)
library(plot.matrix)
library(ggplot2)
library(reshape2)
library(grid)
library(MASS)
library(abind)

set.seed(1234)

# ----------------------------
# Simulation parameters, change targets and noise to the desired simulation scenario 
# ----------------------------
coefs       <- c(rep(1, 10), rep(0, 20)) # targets + noise
nsim        <- 100
nvar        <- 30
simsd       <- 5 * 10
n_subjects  <- 15
targets     <- 10

# Initialize matrices to store results
gllcount           <- matrix(0, nrow = 1, ncol = nvar)
zell_bf_sure_mat   <- matrix(0, nrow = nvar, ncol = nsim)
zell_bf_sqrt_mat   <- matrix(0, nrow = nvar, ncol = nsim)
pos_median_matrix   <- matrix(0, nrow = 6 * nvar, ncol = nsim)
pos_mean_matrix     <- matrix(0, nrow = 6 * nvar, ncol = nsim)
pos_median          <- matrix(0, nrow = 1, ncol = nvar)

# ----------------------------
# Zellner g Prior Function
# ----------------------------
zell_g <- function(Y, X, b0, g = NULL) {
  n <- length(Y)
  p <- ncol(X)
  beta_hat <- solve(t(X) %*% X) %*% (t(X) %*% Y)
  Y_hat <- X %*% beta_hat
  Y0 <- X %*% b0
  sigmasq_hat <- norm(Y - X %*% beta_hat, type = "2")^2 / (n - p)
  
  if (is.null(g)) g <- sqrt(n)
  else if (g == "SURE") g <- (norm(Y_hat - Y0, type = "2")^2 / (p * sigmasq_hat)) - 1
  
  V_star <- (g / (1 + g)) * solve(t(X) %*% X)
  b_star <- (1 / (1 + g)) * b0 + (g / (1 + g)) * beta_hat
  linear_model <- lm(Y ~ 0 + X)
  r_sq <- summary(linear_model)$r.squared
  bf <- (1 + g)^((n - p - 1) / 2) * (1 + g * (1 - r_sq))^(-(n - 1) / 2)
  
  return(list(mean = b_star, sd = V_star * sigmasq_hat, g = g, V_star = V_star, bf = bf))
}

# ----------------------------
# Delta X calculation
# ----------------------------
get_delta_x <- function(x, n, p, t) {
  DXarray <- array(0, dim = c(n, p, t - 1))
  for (k in 1:(t - 1)) {
    DXarray[, , k] <- scale(x[, , k + 1] - x[, , k], center = FALSE) / sqrt(n - 1)
  }
  
  DX <- matrix(0, nrow = n * (t - 1), ncol = sum(1:(t - 1)) * p)
  rowvec <- 1:n
  colvec <- 1:p
  for (i in 1:(t - 1)) {
    for (j in 1:i) {
      DX[rowvec, colvec] <- DXarray[, , j]
      colvec <- colvec + p
    }
    rowvec <- rowvec + n
  }
  
  return(list(DX = DX, DXarray = DXarray))
}

# ----------------------------
# 1) Simulation Loop (No Plots)
# ----------------------------
for (sim in 1:nsim) {
  print(paste("Simulation:", sim))
  
  # Covariance & X matrices
  Sigma <- diag(sqrt(runif(nvar, 5, 5)))
  randnorm <- t(MASS::mvrnorm(n_subjects, mu = runif(nvar, 20, 30), Sigma = Sigma))
  dmu2 <- rep(0, nvar); dmu2[1:targets] <- seq(5, 10, length.out = targets)
  
  x1 <- matrix(randnorm, ncol = nvar)
  x2 <- x1 + matrix(MASS::mvrnorm(n_subjects, mu = dmu2, Sigma = Sigma), ncol = nvar)
  x3 <- x2 + matrix(MASS::mvrnorm(n_subjects, mu = rep(0, nvar), Sigma = Sigma), ncol = nvar)
  x4 <- x3 + matrix(MASS::mvrnorm(n_subjects, mu = rep(0, nvar), Sigma = Sigma), ncol = nvar)
  xpro <- abind::abind(x1, x2, x3, x4, along = 3)
  x <- rbind(x1, x2, x3, x4)
  
  # Generate Y values
  y1 <- abs(rnorm(n_subjects, 100, simsd))
  y2 <- MASS::mvrnorm(1, y1 + colSums(coefs * abs(t(x2 - x1))), diag(rep(simsd, n_subjects)))
  y3 <- MASS::mvrnorm(1, y2 + colSums(coefs * abs(t(x3 - x2))), diag(rep(simsd, n_subjects)))
  y4 <- MASS::mvrnorm(1, y3 + colSums(coefs * abs(t(x4 - x3))), diag(rep(simsd, n_subjects)))
  
  y1 <- y1 + rnorm(n_subjects, 0, 10)
  y2 <- y2 + rnorm(n_subjects, 0, 10)
  y3 <- y3 + rnorm(n_subjects, 0, 10)
  y4 <- y4 + rnorm(n_subjects, 0, 10)
  
  # Consecutive differences for regression
  X <- get_delta_x(xpro, n_subjects, nvar, 4)$DX
  DY2 <- y2 - y1; DY3 <- y3 - y2; DY4 <- y4 - y3
  Y <- c(DY2, DY3, DY4)
  
  # -----------------------------------------------
  # WARNING: The following BGLSS code is commented out
  # because it is extremely slow without parallelization.
  # It calculates posterior medians and coefficient counts.
  # Uncomment at your own risk if you have patience and resources.
  # -----------------------------------------------
  # mod = BGLSS(Y, X[, c(rep(rep((1:nvar), each = 6) + rep(seq(0, (6*nvar - nvar), nvar), nvar))],
  #             niter = 10000, burnin = 5000, group_size = rep(6, nvar),
  #             a = 1, b = 1, num_update = 100, niter.update = 100, verbose = FALSE,
  #             alpha = 0.1, gamma = 0.1, pi_prior = TRUE, pi = 0.5,
  #             update_tau = TRUE, option.weight.group = FALSE,
  #             option.update = "global", lambda2_update = NULL)
  # pos_median_matrix[, sim] <- mod$pos_median
  # for (j in 1:nvar) {
  #   if (sum(abs(mod$coef[(j-1)*6 + 1:6])) > 0) { gllcount[1, j] = gllcount[1, j] + 1 }
  #   if (sum(abs(mod$pos_median[(j-1)*6 + 1:6])) > 0) { pos_median[1, j] = pos_median[1, j] + 1 }
  # }
  
  # Calculate Bayes Factors for each variable
  zell_bf_sure <- zell_bf_sqrt <- rep(0, nvar)
  for (j in 1:nvar) {
    X_uni <- X[, j + nvar * (0:5)]
    bf_sure <- zell_g(Y, X_uni, b0 = rep(0, 6), g = "SURE")
    bf_sqrt <- zell_g(Y, X_uni, b0 = rep(0, 6))
    zell_bf_sure[j] <- as.numeric(bf_sure$bf)
    zell_bf_sqrt[j] <- as.numeric(bf_sqrt$bf)
  }
  zell_bf_sure_mat[, sim] <- zell_bf_sure
  zell_bf_sqrt_mat[, sim] <- zell_bf_sqrt
  
  # Save first simulation trajectory for plotting, can change to different simulations
  if (sim == 1) {
    saveRDS(list(x = x, y = c(y1, y2, y3, y4), id = factor(rep(1:n_subjects, 4))),
            file = "trajectory_data.rds")
  }
}

# Save results, make sure to change csv file to reflect simulation scenario
write.csv(zell_bf_sure_mat, "zell_bf_sure_mat_30.csv", row.names = FALSE)
write.csv(zell_bf_sqrt_mat, "zell_bf_sqrt_mat_30.csv", row.names = FALSE)

# ----------------------------
# 2) Plot Trajectories (Separate)
# ----------------------------
traj_data <- readRDS("trajectory_data.rds")
x <- traj_data$x
y <- traj_data$y
id <- traj_data$id
meltY <- data.frame(y = y, X = rep(c(1,3,5,15), each = n_subjects), id = id)

# Example metabolite plot
j <- 25
meltY$x <- x[, j]
ggplot(data = meltY, aes(x = as.factor(X), y = x, group = id)) +
  geom_line(col = '#00BFC4') +
  xlab('Days') + ylab('Log 2 Abundance') +
  ggtitle("Noise Metabolite") +
  scale_x_discrete(breaks = c(1,3,5,15), labels = c('0','2','4','14')) +
  theme(plot.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))

# Overall trajectory plot
meltY$y <- y
ggplot(data = meltY, aes(x = as.factor(X), y = y, group = id)) +
  geom_line(col = '#00BFC4') +
  xlab('Days') + ylab('TTP') +
  ggtitle("Mycobacterial Load Trajectories") +
  scale_x_discrete(breaks = c(1,3,5,15), labels = c('0','2','4','14')) +
  theme(plot.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))
