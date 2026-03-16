# leave-one-out cross-validation to select bandwidth
library(KernSmooth)

w.fn <- function(x.vals, x, bw) {
  w.avals <- NULL
  for (x.val in x.vals) {
    x.std <- (x - x.val) / bw
    kern.std <- dnorm(x.std) / bw
    w.avals <- c(w.avals, mean(x.std ^ 2 * kern.std) * (dnorm(0) / bw) /
                   (mean(kern.std) * mean(x.std ^ 2 * kern.std) - 
                      mean(x.std * kern.std) ^ 2))
  }
  return(w.avals/n)
}

hatvals <- function(x.vals, bw) {
  approx(x.vals, w.fn(x.vals, x, bw), xout=x)$y
}

cts.eff <- function(x, out, bw) {
  approx(locpoly(x = x, y = out, bandwidth = bw), xout=x)$y
}

h.fn <- function(h, x.vals, x, out) {
  hats <- hatvals(x.vals, h)
  mean(((out - cts.eff(x, out, bw=h)) / (1 - hats)) ^ 2)
}

approx.fn <- function(x, y, z) predict(smooth.spline(x, y), x = z)$y

decile_heterogeneity <- function(var, res, xlab) {
  dt <- as.data.table(res)
  dt[, dec := .bincode(get(var), 
                        quantile(get(var), seq(0, 1, 0.1)),
      include.lowest = TRUE)]
  z_te <- dt[, list(zte = mean(cate), dev = sd(cate) / sqrt(.N), cnt = .N), 
              by = c("dec")]
  ggplot(z_te, aes(x = dec, y = zte)) +
    geom_point() + theme_bw() + geom_line() +
    geom_ribbon(aes(ymin=zte - 1.96 * dev, ymax = zte + 1.96 * dev), alpha = 0.3) +   
    labs(x = xlab,
         y = "z-TE")
}
