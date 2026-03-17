library(mice)
library(qgcomp)
source('cleaning.R')

set.seed(123)

sum(df$cvd == 0, na.rm = TRUE)

df$cvd <- ifelse(df$cvd == 0, 0.1 / df$population, df$cvd)

N <- nrow(df)

lod_Y <- ifelse(is.na(df$cvd), 
                     9 / df$population, 
                     NA)

lodlist <- lapply(df, function(x) rep(NA, N))
lodlist$cvd <- lod_Y

init <- mice(df, maxit = 0)
meth <- init$method
meth["cvd"] <- "tobit"

mi_dat <- mice(df, 
               method = meth, 
               lod = lodlist, 
               m = 5,            
               maxit = 10,       
               seed = 123)

complete_data_1 <- complete(mi_dat, 1)

max_deaths_imputed <- max(complete_data_1$cvd[is.na(df$cvd)] * 
                            df$population[is.na(df$cvd)])
print(max_deaths_imputed)                  
                   
