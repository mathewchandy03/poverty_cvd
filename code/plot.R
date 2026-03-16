source('variation_analysis.R')

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

