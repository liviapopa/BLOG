library(MBSGS)
library(MASS)
library(abind)
library(doParallel)
library(foreach)

cat("Started.\n")
set.seed(1234)

### ==== Simulation constants ====
#change this for however many variables & targets you are simulating
n_columns = 350
targets = 50

#c(targets,noise)
coefs = c(rep(1,targets), rep(0,n_columns-targets))
coefs3 = coefs4 = c(rep(0,targets), rep(0,n_columns-targets))
simsd = 5 * targets

#number simulations
nsim = 100
gllcount = llcount = matrix(0, nrow = 1, ncol = n_columns)
zell_bf_sqrt_mat = matrix(0, nrow = n_columns, ncol = nsim)
zell_bf_sure_mat = matrix(0, nrow = n_columns, ncol = nsim)
pos_median_matrix = matrix(0, nrow = n_columns*6, ncol = nsim)
pos_mean_matrix = matrix(0, nrow = n_columns*6, ncol = nsim)
pos_median = matrix(0, nrow = 1, ncol = n_columns)

### ==== Functions ====

# Zellner's g-prior (returns posterior mean, V*, Bayes factor etc.)
zell_g <- function(Y, X, b0, g = NULL) {
  n <- length(Y)
  p <- ncol(X)
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
  linear_model <- lm(Y ~ 0 + X)
  r_sq <- summary(linear_model)$r.squared
  bf <- (1 + g)^((n - p - 1) / 2) * (1 + g * (1 - r_sq))^(-(n - 1) / 2)
  
  return(list(mean = b_star, sd = V_star * sigmasq_hat, g = g, md = md, V_star = V_star, bf = bf))
}

# Build delta-X matrix and 3D array for t time points
get_delta_x <- function(x, n, p, t) {
  DXarray <- array(0, dim = c(n, p, t - 1))
  for (k in 1:(t - 1)) {
    DXarray[, , k] <- scale(x[, , k + 1] - x[, , k], center = F) / sqrt(n - 1)
  }
  colnames(DXarray) <- colnames(x)
  
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

### ==== Parallel setup ====
cl <- makeCluster(detectCores() - 2) # or any other number of cores you want to use
registerDoParallel(cl)

results <- foreach(sim = 1:nsim, .packages = c("MBSGS", "MASS", "abind"), 
                   .export = c("zell_g", "get_delta_x", "coefs", "coefs3", "coefs4", "simsd"), 
                   .combine = 'c') %dopar% { #results[[1]]$gll_count
                     
                     # Simulation parameters change this to however many targets and noise you want
                     n_subjects = 15
                     n_col = 350
                     targets = 50
                     
                     # Generate random data for simulation
                     Sigma = diag(sqrt(runif(n_col, 5, 5)))
                     randnorm = t(MASS::mvrnorm(n_subjects, mu = runif(n_col, 20, 30), Sigma = Sigma))
                     
                     dmu2 = dmu3 = dmu4 = rep(0, n_col)
                     dmu2[1:targets] = seq(5, 10, length.out = targets)
                     
                     x1 = matrix(randnorm, ncol = n_col)
                     x2 = x1 + matrix(MASS::mvrnorm(n_subjects, mu = dmu2, Sigma = Sigma), ncol = n_col)
                     x3 = x2 + matrix(MASS::mvrnorm(n_subjects, mu = dmu3, Sigma = Sigma), ncol = n_col)
                     x4 = x3 + matrix(MASS::mvrnorm(n_subjects, mu = dmu4, Sigma = Sigma), ncol = n_col)
                     
                     xpro = abind::abind(x1, x2, x3, x4, along = 3)
                     x = rbind(x1, x2, x3, x4)
                     
                     y1 = abs(rnorm(n_subjects, 100, simsd))
                     y2 = MASS::mvrnorm(1, y1 + colSums(coefs * abs(t(x2 - x1))), diag(rep(simsd, n_subjects)))
                     y3 = MASS::mvrnorm(1, y2 + colSums(coefs * abs(t(x2 - x1))) + colSums(coefs3 * abs(t(x3 - x2))), diag(rep(simsd, n_subjects)))
                     y4 = MASS::mvrnorm(1, y3 + colSums(coefs * abs(t(x2 - x1))) + colSums(coefs3 * abs(t(x3 - x2))) + colSums(coefs4 * abs(t(x4 - x3))), diag(rep(simsd, n_subjects)))
                     
                     time = factor(rep(4:1, each = n_subjects))
                     id = factor(rep(1:n_subjects, 4))
                     ypro = (cbind(y1, y2, y3, y4))
                     y = (c(y1, y2, y3, y4))
                     
                     X = get_delta_x(xpro, n_subjects, n_col, 4)$DX
                     
                     DY4 = y4 - y3
                     DY3 = y3 - y2
                     DY2 = y2 - y1
                     Y = c(DY2, DY3, DY4)
                     
                     X_consecutive = X[, c(rep(rep((1:n_col), each = 6) + rep(seq(0, (6*n_col - n_col), n_col), n_col)))]
                     
                     # Run BGLSS model
                     mod = BGLSS(Y, X_consecutive, niter = 10000, burnin = 5000, group_size = rep(6, n_col), a = 1, 
                                 b = 1, num_update = 100, niter.update = 100, verbose = FALSE,
                                 alpha = 0.1, gamma = 0.1, pi_prior = FALSE, pi = 0.99, 
                                 update_tau = TRUE, option.weight.group = FALSE,
                                 option.update = "global", lambda2_update = NULL)
                     
                     # Initialize local variables to store results for this simulation
                     sim_gllcount <- matrix(0, nrow = 1, ncol = n_col)
                     sim_pos_median <- matrix(0, nrow = 1, ncol = n_col)
                     sim_zell_bf_sqrt <- rep(0, n_col)
                     sim_zell_bf_sure <- rep(0, n_col)
                     sim_pos_median_matrix <- rep(0,n_col)
                     sim_pos_mean_matrix <- rep(0,n_col)
                     
                     sim_pos_median_matrix <- mod$pos_median 
                     sim_pos_mean_matrix <- mod$pos_mean
                     # Check if coefficients are non-zero, update counters for gllcount and pos_median
                     for (j in 1:n_col) {
                       if (sum(abs(mod$coef[(j - 1) * 6 + 1:6])) > 0) {
                         sim_gllcount[1, j] <- sim_gllcount[1, j] + 1
                       }
                       if (sum(abs(mod$pos_median[(j - 1) * 6 + 1:6])) > 0) {
                         sim_pos_median[1, j] <- sim_pos_median[1, j] + 1
                       }
                     }
                     
                     # Calculate Bayes factors for SURE and regular methods
                     for (j in 1:n_col) {
                       X_uni <- X[, j + n_col * (0:5)]
                       bf_sure <- zell_g(Y, X_uni, b0 = rep(0, 6), g = "SURE")
                       bf_sqrt <- zell_g(Y, X_uni, b0 = rep(0, 6))
                       sim_zell_bf_sure[j] <- as.numeric(bf_sure$bf)
                       sim_zell_bf_sqrt[j] <- as.numeric(bf_sqrt$bf)
                     }
                     
                     # Return results for this simulation
                     list(gllcount = sim_gllcount, pos_median = sim_pos_median, 
                          zell_bf_sqrt = t(sim_zell_bf_sqrt), zell_bf_sure = t(sim_zell_bf_sure),
                          pos_median_matrix = t(sim_pos_median_matrix), pos_mean_matrix = t(sim_pos_mean_matrix))
                   }#foreach

# Stop the parallel cluster
stopCluster(cl)

# Aggregate results outside the loop
gllcount_new = (matrix(NA, nrow = n_columns, ncol = 1))
zell_bf_sqrt_mat_new = (matrix(NA, nrow = n_columns, ncol = 1))
zell_bf_sure_mat_new = (matrix(NA, nrow = n_columns, ncol = 1))
pos_median_matrix_new = (matrix(NA, nrow = n_columns*6, ncol = 1))
pos_mean_matrix_new = (matrix(NA, nrow = n_columns*6, ncol = 1))
pos_median_new = (matrix(NA, nrow = n_columns, ncol = 1))

for (i in 1:length(results)) {
  cat("i = ",i,"\n")
  results_i = t(results[[i]])
  cat("dim results i: ",dim(results_i), "\n")
  if(i%%6==1){
    gllcount_new = cbind(gllcount_new, results_i)
  }
  if(i%%6==2){
    pos_median_new = cbind(pos_median_new, results_i)
  }
  if(i%%6==3)
    zell_bf_sqrt_mat_new = cbind(zell_bf_sqrt_mat_new, results_i)
  if(i%%6==4)
    zell_bf_sure_mat_new = cbind(zell_bf_sure_mat_new, results_i)
  if(i%%6==5)
    pos_median_matrix_new = cbind(pos_median_matrix_new, results_i)
  if(i%%6==0)
    pos_mean_matrix_new = cbind(pos_mean_matrix_new, results_i)
  
}#for

# Finalizing the output
gllcount_new = rowSums(gllcount_new, na.rm = TRUE)
pos_median_new = rowSums(pos_median_new,na.rm = TRUE)
zell_bf_sqrt_mat_new = zell_bf_sqrt_mat_new[,-1] #remove NA first col
zell_bf_sure_mat_new = zell_bf_sure_mat_new[,-1]
pos_median_matrix_new = pos_median_matrix_new[,-1]
pos_mean_matrix_new = pos_mean_matrix_new[,-1]

# Write the final results to CSV
#change to whatever signal to noise ratio you want to save
write.csv(gllcount_new, "gllcount_350_pi99.csv", row.names = FALSE)
write.csv(pos_median_matrix_new, "pos_median_mat_350_pi99.csv", row.names = FALSE)
write.csv(zell_bf_sqrt_mat_new, "zell_bf_sqrt_mat_350.csv", row.names = FALSE)
write.csv(zell_bf_sure_mat_new, "zell_bf_sure_mat_350.csv", row.names = FALSE)
write.csv(pos_mean_matrix_new, "pos_mean_matrix_350_pi99.csv", row.names = FALSE)
write.csv(pos_median_new, "pos_median350_pi99.csv", row.names = FALSE)

gll_count350 <- read.csv("gllcount_350_pi99.csv")
View(gll_count350)

pos_median100 <- read.csv("pos_median350_pi99.csv")
View(pos_median100)
