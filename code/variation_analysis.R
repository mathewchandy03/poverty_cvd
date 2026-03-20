library(SuperLearner)
library(dplyr)
library(tidyr)
library(ggplot2)
library(KernSmooth)
source('imputation.R')
source('functions.R')

set.seed(123)

m_imputations <- 5
n <- nrow(df)

tv_matrix <- matrix(NA, nrow = 100, ncol = m_imputations)
te_matrix <- matrix(NA, nrow = 100, ncol = m_imputations)
se_matrix <- matrix(NA, nrow = 100, ncol = m_imputations)
# New matrix to store the standard errors of the Total Effect
te_se_matrix <- matrix(NA, nrow = 100, ncol = m_imputations) 
tv_se_matrix <- matrix(NA, nrow = 100, ncol = m_imputations) 

# Extract covariates/treatment/response
l <- as.matrix(df %>% select(nonwhite_rate, age65plus, rurality, medicaid))
x <- df$poverty_rate
y <- df$cvd

# set up evaluation points & matrices for predictions
n <- nrow(l)
x.min <- min(x)
x.max <- max(x)
# x.min <- quantile(x, 0.05)
# x.max <- quantile(x, 0.95)
x.vals <- seq(x.min, x.max, length.out=100)

# Define kernel function for the standard error calculation
kern <- dnorm 

for (i in 1:m_imputations) {
  current_data <- complete(mi_dat, i)
  
  y <- current_data$cvd
  x <- current_data$poverty_rate
  l <- as.matrix(current_data %>% 
                   select(nonwhite_rate, age65plus, rurality, medicaid))
  w <- current_data$population
  
  lx.new <- rbind(cbind(l, x), cbind(l[rep(1:n, length(x.vals)), ], 
                                     x = rep(x.vals, each = n)))
  l.new <- lx.new[,-dim(lx.new)[2]]
  
  # convert l to data.frame to suppress the X is not a data frame warning
  l <- as.data.frame(l) 
  lx.new <- as.data.frame(lx.new)
  l.new <- as.data.frame(l.new)
  
  # fit super learner
  sl.lib <- c("SL.xgboost")
  pimod <- SuperLearner(Y = x, X = l, SL.library = sl.lib, newX = l.new, 
                        obsWeights = w)
  pimod.vals <- pimod$SL.predict
  sq.res <- (x - pimod.vals[1:n]) ^ 2
  pi2mod <- SuperLearner(Y = sq.res, X = l, SL.library = sl.lib, newX = l.new,
                         obsWeights = w)
  pi2mod.vals <- pmax(pi2mod$SL.predict, 1e-4)
  mumod <- SuperLearner(Y=y, X=cbind(l,x), SL.library=sl.lib,
                        newX=lx.new, family=gaussian, obsWeights = w)
  muhat.vals <- mumod$SL.predict
  
  # construct estimated p/varpi and mu/m values
  x.std <- (lx.new$x - pimod.vals) / sqrt(pi2mod.vals)
  pihat.vals <- approx.fn(density(x.std[1:n])$x, density(x.std[1:n])$y,
                          x.std)
  pihat <- pihat.vals[1:n]
  pihat.mat <- matrix(pihat.vals[-(1:n)], nrow = n, ncol = length(x.vals))
  varpihat <- approx.fn(x.vals, apply(pihat.mat, 2, mean), x)
  varpihat.mat <- matrix(rep(apply(pihat.mat, 2, mean), n), byrow = T, nrow = n)
  
  muhat <- muhat.vals[1:n]
  muhat.mat <- matrix(muhat.vals[-(1:n)], nrow = n, ncol = length(x.vals))
  mhat <- approx.fn(x.vals, apply(muhat.mat, 2, mean), x)
  mhat.mat <- matrix(rep(apply(muhat.mat, 2, mean), n), byrow = T, nrow = n)
  
  dens_ratio <- pihat / varpihat
  
  # dens_ratio <- pmax(dens_ratio, 0.05)
  
  # form adjusted/pseudo outcome xi
  pseudo.out <- (y - muhat) / dens_ratio + mhat
  
  # choose bandwidth for total variation
  h.opt_tv <- optimize(h.fn, c(0.01, 0.5), tol=0.01, 
                       x.vals = x.vals, x = x, out = y)$minimum
  tv_matrix[, i] <- approx(locpoly(x, y, bandwidth=h.opt_tv), xout = x.vals)$y
  
  # choose bandwidth for total effect
  h.opt_te <- optimize(h.fn, c(0.01, 0.5), tol=0.01, 
                       x.vals = x.vals, x = x, out = pseudo.out)$minimum
  te_matrix[, i] <- approx(locpoly(x, pseudo.out, bandwidth=h.opt_te), xout = x.vals)$y
  
  se_matrix[, i] <- tv_matrix[, i] - te_matrix[, i]
  
  # estimate sandwich-style pointwise confidence band
  tv_se_imputation <- NULL
  te_se_imputation <- NULL
  
  for (x.val in x.vals) {
    x.std.tv <- (x - x.val) / h.opt_tv
    kern.std.tv <- (kern(x.std.tv) / h.opt_tv) / h.opt_tv
    
    beta.tv <- coef(lm(y ~ x.std.tv, weights = w * kern.std.tv))
    
    Dh.tv <- matrix(c(mean(w * kern.std.tv), mean(w * kern.std.tv * x.std.tv),
                      mean(w * kern.std.tv * x.std.tv), mean(w * kern.std.tv * x.std.tv^2)), nrow=2)
    
    iffn.tv <- w * kern.std.tv * (y - beta.tv[1] - beta.tv[2] * x.std.tv)
    sigma.tv <- cov(t(solve(Dh.tv) %*% rbind(iffn.tv, x.std.tv * iffn.tv)))
    
    tv_se_imputation <- c(tv_se_imputation, sqrt(sigma.tv[1,1] / n))
    
    x.std.val <- (x - x.val) / h.opt_te
    kern.std <- (kern(x.std.val) / h.opt_te) / h.opt_te
    
    beta <- coef(lm(pseudo.out ~ x.std.val, weights = w * kern.std))
    
    Dh <- matrix(c(mean(w * kern.std), mean(w * kern.std * x.std.val),
                   mean(w * kern.std * x.std.val), mean(w * kern.std * x.std.val^2)), nrow=2)
    
    kern.mat <- matrix(rep(kern((x.vals - x.val) / h.opt_te) / h.opt_te, n), byrow=T, nrow=n)
    g2 <- matrix(rep((x.vals - x.val) / h.opt_te, n), byrow=T, nrow=n)
    
    intfn1.mat <- kern.mat * (muhat.mat - mhat.mat) * varpihat.mat
    intfn2.mat <- g2 * kern.mat * (muhat.mat - mhat.mat) * varpihat.mat
    
    grid_diff <- (x.vals[-1] - x.vals[-length(x.vals)]) / 2
    int1 <- apply(matrix(rep(grid_diff, n), byrow=T, nrow=n) * intfn1.mat[, -1] + intfn1.mat[, -length(x.vals)], 1, sum)
    int2 <- apply(matrix(rep(grid_diff, n), byrow=T, nrow=n) * intfn2.mat[, -1] + intfn2.mat[, -length(x.vals)], 1, sum)
    
    iffn.te.1 <- w * kern.std * (pseudo.out - beta[1] - beta[2] * x.std.val) + int1
    iffn.te.2 <- x.std.val * w * kern.std * (pseudo.out - beta[1] - beta[2] * x.std.val) + int2
    
    sigma <- cov(t(solve(Dh) %*% rbind(iffn.te.1, iffn.te.2)))
    
    te_se_imputation <- c(te_se_imputation, sqrt(sigma[1,1] / n))
  }
  te_se_matrix[, i] <- te_se_imputation
  tv_se_matrix[, i] <- tv_se_imputation
}

# --- Aggregate using Rubin's Rules ---
final_tv_curve <- rowMeans(tv_matrix, na.rm = TRUE)
final_te_curve <- rowMeans(te_matrix, na.rm = TRUE)
final_se_curve <- rowMeans(se_matrix, na.rm = TRUE)

W_var_te <- rowMeans(te_se_matrix^2, na.rm = TRUE)
B_var_te <- apply(te_matrix, 1, var, na.rm = TRUE)
total_te_se <- sqrt(W_var_te + (1 + 1/m_imputations) * B_var_te)

W_var_tv <- rowMeans(tv_se_matrix^2, na.rm = TRUE)
B_var_tv <- apply(tv_matrix, 1, var, na.rm = TRUE)
total_tv_se <- sqrt(W_var_tv + (1 + 1/m_imputations) * B_var_tv)

df_ribbons <- data.frame(
  x = x.vals * 100,
  te_ci_lower = final_te_curve - 1.96 * total_te_se,
  te_ci_upper = final_te_curve + 1.96 * total_te_se,
  tv_ci_lower = final_tv_curve - 1.96 * total_tv_se,
  tv_ci_upper = final_tv_curve + 1.96 * total_tv_se
)

df_lines <- data.frame(
  x = x.vals * 100,
  tv = final_tv_curve,
  te = final_te_curve,
  se = final_se_curve
) %>% 
  pivot_longer(cols = c(tv, te, se), names_to = "Quantity", values_to = "value")



ggplot() +
  geom_ribbon(data = df_ribbons, aes(x = x, ymin = te_ci_lower, ymax = te_ci_upper), 
              alpha = 0.2, fill = "#F8766D") + 
  geom_ribbon(data = df_ribbons, aes(x = x, ymin = tv_ci_lower, ymax = tv_ci_upper), 
              alpha = 0.2, fill = "#00BFC4") + 
  geom_line(data = df_lines %>% filter(Quantity != 'se'), aes(x = x, y = value, color = Quantity), linewidth = 1) +
  labs(x = "Poverty Rate (%)",
       y = "Age-adjusted Cardiovascular Disease \n Mortality per 100,000") +
  theme_minimal() +
  scale_color_discrete(labels=c("Total Effect", "Total Variation"))

ggsave("../manuscript/figures/variation_analysis.pdf", width = 6, height = 4)

