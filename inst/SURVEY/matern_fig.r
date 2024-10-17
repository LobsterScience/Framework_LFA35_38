#visualizing the matern

range <- 114 #< try changing this

kappa <- sqrt(8) / range
distance <- seq(0.001, 300, length.out = 100)
correlation <- kappa * distance * besselK(kappa * distance, nu = 1)
plot(distance, correlation, type = "l", ylim = c(0, 1))
abline(h = 0.13, lty = 2) # 0.13 correlation
abline(v = range, lty = 2)