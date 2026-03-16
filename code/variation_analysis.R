library(SuperLearner)
source('imputation.R')
source('functions.R')

m_imputations <- 5
n <- nrow(df)

tv_matrix <- matrix(NA, nrow = 100, ncol = m_imputations)
te_matrix <- matrix(NA, nrow = 100, ncol = m_imputations)
se_matrix <- matrix(NA, nrow = 100, ncol = m_imputations)

# Extract covariates/treatment/response
l <- as.matrix(df %>% select(nonwhite_rate, age65plus, rurality, medicaid))
x <- df$poverty_rate
y <- df$cvd

# set up evaluation points & matrices for predictions
n <- nrow(l)
x.min <- min(x)
x.max <- max(x)
x.vals <- seq(x.min, x.max, length.out=100)

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
  # varpihat.mat <- matrix(rep(apply(pihat.mat, 2, mean), n), byrow = T, nrow = n)
  muhat <- muhat.vals[1:n]
  muhat.mat <- matrix(muhat.vals[-(1:n)], nrow = n, ncol = length(x.vals))
  mhat <- approx.fn(x.vals, apply(muhat.mat, 2, mean), x)
  # mhat.mat <- matrix(rep(apply(muhat.mat, 2, mean), n), byrow = T, nrow = n)
  
  # form adjusted/pseudo outcome xi
  pseudo.out <- (y - muhat) / (pihat / varpihat) + mhat
  
  # choose bandwidth for total variation
  
  h.opt <- optimize(h.fn, c(0.01, 0.5), tol=0.01, 
                    x.vals = x.vals, x = x, out = y)$minimum
  
  # estimate variation curve with optimal bandwidth
  
  tv_matrix[, i] <- approx(locpoly(x, y, bandwidth=h.opt), xout = x.vals)$y
  
  # choose bandwidth for total effect
  
  h.opt <- optimize(h.fn, c(0.01, 0.5), tol=0.01, 
                    x.vals = x.vals, x = x, out = pseudo.out)$minimum
  
  # estimate effect curve with optimal bandwidth
  
  te_matrix[, i] <- approx(locpoly(x, pseudo.out, bandwidth=h.opt), xout = x.vals)$y
  
  # # choose bandwidth for spurious effect
  # 
  # h.opt <- optimize(h.fn, c(0.01, 0.5), tol=0.01, 
  #                   x.vals = x.vals, x = x, out = pseudo.out - y)$minimum
  # 
  # # estimate effect curve with optimal bandwidth
  # 
  # se_matrix[, i] <- approx(locpoly(x, pseudo.out - y, bandwidth=h.opt), xout = x.vals)$y
  se_matrix[, i] <- tv_matrix[, i] - te_matrix[, i]
}

final_tv_curve <- rowMeans(tv_matrix, na.rm = TRUE)
final_te_curve <- rowMeans(te_matrix, na.rm = TRUE)
final_se_curve <- rowMeans(se_matrix, na.rm = TRUE)

df_wide <- data.frame(
  x = x.vals * 100,
  tv = final_tv_curve ,
  te = final_te_curve,
  se = final_se_curve
)

df_long <- pivot_longer(df_wide, cols = c(tv, te, se), 
                        names_to = "Quantity", values_to = "value")

ggplot(df_long, aes(x = x, y = value, color = Quantity)) +
  geom_line() +
  labs(title = "The Effect of Poverty Rate on Age-adjusted Cardiovascular 
       Disease Mortality per 1000",
       x = "Poverty Rate (%)",
       y = "Age-adjusted Cardiovascular Disease 
       Mortality per 1000") +
  theme_minimal() +
  scale_color_discrete(labels=c("Spurious Effect", "Total Effect", 
                                "Total Variance"))




