source('cleaning.R')

set.seed(123)

n <- nrow(df)
X <- df$poverty_rate
C <- df$nonwhite_rate
A <- df$age65plus
P <- df$population
M <- df$uninsurance
R <- df$rurality
G <- df$medicaid

w1 <- sum(A * P) / sum(P)
w2 <- sum((1 - A) * P) / sum(P)

beta <- 15

nrep <- 100
res <- NULL
for (i in seq_len(nrep)) {
  print(i)
  Y1 <- rnorm(n, mean = (1.0 + 15 * X + 5 * C + 4 * M - 0.2 * R - 0.2 * G), sd = 2) 
  Y1 <- pmax(0, Y1) 
  Y2 <- rnorm(n, mean = (10.0 + 15 * X + 5 * C + 4 * M - 0.2 * R - 0.2 * G), sd = 4) 
  Y2 <- pmax(0, Y2)
  Y_raw <- ( (P * (1 - A) * Y1 / 1000) + (P * A * Y2 / 1000) ) / P * 1000 
  # unadjusted CVD mortality
  Y <- sum((1 - A) * P) / sum(P) * Y1 + sum(A * P) / sum(P) * Y2 # age-adjusted
  Ry <- (P * Y_raw / 1000) < 10
  
  # generate the available data P(V*, R)
  Ys <- Y
  Ys[as.logical(Ry)] <- NA
  data <- data.frame(X, C, A, P, M, R, G, Ys)
  
  # apply MICE with 5 imputations
  mi_dat <- mice(data = data, m = 5, print=FALSE)
  mirh <- c()
  for (mi in seq_len(mi_dat$m)) {
    current_data <- complete(mi_dat, action = mi)
    fit <- lm(Ys ~ X + C + M + R + G, data = current_data)
    mirh <- c(mirh, coef(fit)["X"])
  }
  
  # apply left censored MICE with 5 imputations
  data$Ys <- ifelse(data$Ys == 0, 100 / data$P, data$Ys)
  lod_Y <- ifelse(is.na(data$Ys), 
                  9 / data$P * 1000, 
                  NA)
  lodlist <- lapply(data, function(x) rep(NA, n))
  lodlist$Ys <- lod_Y
  init <- mice(data, maxit = 0)
  meth <- init$method
  meth["Ys"] <- "tobit"  
  mi_dat <- mice(data, 
                 method = meth, 
                 lod = lodlist, 
                 m = 5,            
                 maxit = 10,       
                 seed = i, 
                 print=FALSE)
  mirh_lc <- c()
  for (mi in seq_len(mi_dat$m)) {
    current_data <- complete(mi_dat, action = mi)
    fit <- lm(Ys ~ X + C + M + R + G, data = current_data)
    mirh_lc <- c(mirh_lc, coef(fit)["X"])
  }
  
  beta_mice <- mean(mirh)
  beta_mice_lc <- mean(mirh_lc)
  complete_data <- data[complete.cases(data), ]
  fit <- lm(Ys ~ X + C + M + R + G, data = complete_data)
  beta_cca <- coef(fit)["X"]
  res <- rbind(
    res, 
    data.frame(mice_bias = beta_mice - beta, mice_lc_bias = beta_mice_lc - beta,
               cca_bias = beta_cca - beta)
  )
}

res <- as.data.table(res)
res_long <- melt(res, measure.vars = 1:3)

# plot the bias distribution over repetitions
ggplot(res_long, aes(x=value, fill = variable)) +
  geom_density(alpha = 0.4) + geom_vline(xintercept = 0, color = "red",
                                         linetype = "dashed") +
  theme_bw() + scale_fill_discrete(name = "Method", 
                                   labels = c("Standard MICE", "Left Censored MICE", "CCA")) +
  xlab("Bias") + ylab("Distribution Density")

ggsave("../manuscript/figures/sensitivity_analysis.pdf", width = 6, height = 4)
