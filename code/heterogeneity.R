library(grf)
library(xtable)
source("imputation.R")
source("functions.R")


m_imputations <- 5
n <- nrow(df)

cate <- matrix(NA, nrow = n, ncol = m_imputations)

for (i in 1:m_imputations) {
  res <- complete(mi_dat, i)
  
  crf <- causal_forest(
    X = as.matrix(res[, c("nonwhite_rate", "age65plus", "rurality", "medicaid")]), # all covariates -- pre-treatment!
    Y = as.numeric(res$cvd), # outcome variable
    W = res$poverty_rate # treatment variable
  )
  
  # get conditional treatment effects (z-TE) for each row
  cate[, i] <- crf$predictions
}

df$cate <- rowMeans(cate, na.rm = TRUE)

print(mean(df$cate))

df$poverty_level <- ifelse(df$poverty_rate > median(df$poverty_rate), 
                            "High Poverty (High Dose)", 
                            "Low Poverty (Low Dose)")

ame_t <- df %>%
  group_by(poverty_level) %>%
  summarize(average_marginal_effect = mean(cate / 100)) 

colnames(ame_t) <- c("Poverty Level", "Average Marginal Effect")

xt <- xtable(ame_t, type = "latex", label = 'tab:ett',
             caption = "Table of average marginal effect for binarized version 
             of poverty level, where poverty levels greater than the median are 
             considered high, and poverty levels lesser than the median 
             are considered low.")

print(xt,
      file = "../manuscript/tables/ett.tex",
      caption.placement = "top",
      include.rownames = FALSE)


decile_heterogeneity("nonwhite_rate", df, 
                     xlab = "Proportion of Non-Whites Decile")

ggsave("../manuscript/figures/z_te_nonwhite_rate.pdf", width = 6, height = 4)

decile_heterogeneity("age65plus", df, 
                     xlab = "Proportion of Age 65+ Decile")

ggsave("../manuscript/figures/z_te_age65plus.pdf", width = 6, height = 4)

dt <- as.data.table(df)
z_te <- dt[, list(zte = mean(cate), dev = sd(cate) / sqrt(.N), cnt = .N), 
           by = c("rurality")]
ggplot(z_te, aes(x = rurality, y = zte)) +
  geom_point() + theme_bw() + geom_line() +
  geom_ribbon(aes(ymin=zte - 1.96 * dev, ymax = zte + 1.96 * dev), alpha = 0.3) +   
  labs(x = "Urban-Rural Classification Scheme",
       y = "z-TE")

ggsave("../manuscript/figures/z_te_rurality.pdf", width = 6, height = 4)

z_te <- dt[, list(zte = mean(cate), dev = sd(cate) / sqrt(.N), cnt = .N), 
           by = c("medicaid")]
ggplot(z_te, aes(x = medicaid, y = zte)) +
  geom_point() + theme_bw() + geom_line() +
  geom_ribbon(aes(ymin=zte - 1.96 * dev, ymax = zte + 1.96 * dev), alpha = 0.3) +   
  labs(x = "Medicaid Expansion Status",
       y = "z-TE")

ggsave("../manuscript/figures/z_te_medicaid.pdf", width = 6, height = 4)







