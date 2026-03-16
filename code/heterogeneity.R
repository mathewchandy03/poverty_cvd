library(grf)
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
df$effect_group <- ifelse(df$cate > median(df$cate), 
                           "High Impact (Steep Slope)", 
                           "Low Impact (Flat Slope)")

df$poverty_level <- ifelse(df$poverty_rate > median(df$poverty_rate), 
                            "High Poverty (High Dose)", 
                            "Low Poverty (Low Dose)")

ame_t <- df %>%
  group_by(poverty_level) %>%
  summarize(average_marginal_effect = mean(cate))

print(ame_t)


decile_heterogeneity("nonwhite_rate", df, 
                     xlab = "Proportion of Non-Whites Decile")

ggsave("../manuscript/figures/z_te_nonwhite_rate.pdf")

decile_heterogeneity("age65plus", df, 
                     xlab = "Proportion of Age 65+ Decile")

ggsave("../manuscript/figures/z_te_age65plus.pdf")

dt <- as.data.table(df)
z_te <- dt[, list(zte = mean(cate), dev = sd(cate) / sqrt(.N), cnt = .N), 
           by = c("rurality")]
ggplot(z_te, aes(x = rurality, y = zte)) +
  geom_point() + theme_bw() + geom_line() +
  geom_ribbon(aes(ymin=zte - 1.96 * dev, ymax = zte + 1.96 * dev), alpha = 0.3) +   
  labs(x = "Urban-Rural Classification Scheme",
       y = "z-TE")

ggsave("../manuscript/figures/z_te_rurality.pdf")

z_te <- dt[, list(zte = mean(cate), dev = sd(cate) / sqrt(.N), cnt = .N), 
           by = c("medicaid")]
ggplot(z_te, aes(x = medicaid, y = zte)) +
  geom_point() + theme_bw() + geom_line() +
  geom_ribbon(aes(ymin=zte - 1.96 * dev, ymax = zte + 1.96 * dev), alpha = 0.3) +   
  labs(x = "Medicaid Expansion Status",
       y = "z-TE")

ggsave("../manuscript/figures/z_te_medicaid.pdf")


ggplot(df, aes(x = nonwhite_rate, color = poverty_level)) +  
  stat_ecdf(linewidth = 1) + 
  theme_bw() +   
  labs(title = "Age Distribution by Poverty Exposure",
       x = "Percent Age 65+",
       y = "Cumulative Proportion")

ggplot(df, aes(x = nonwhite_rate * 100, color = effect_group)) +  
  stat_ecdf(linewidth = 1) + 
  theme_bw() +   
  labs(title = "Non-White Percentage Distribution by Poverty-Vulnerability",
       x = "Non-White Percentage",
       y = "Cumulative Proportion",
       color = "Treatment Effect Size")

ggplot(df, aes(x = age65plus * 100, color = effect_group)) +  
  stat_ecdf(linewidth = 1) + 
  theme_bw() +   
  labs(title = "Age 65+ Percentage Distribution by Poverty-Vulnerability",
       x = "Age 65+ Percentage",
       y = "Cumulative Proportion",
       color = "Treatment Effect Size")




