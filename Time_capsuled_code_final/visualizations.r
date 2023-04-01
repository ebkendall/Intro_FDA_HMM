# First simulate the time series according to some state sequence
set.seed(140)
x = seq(-3,3,by=0.01)

y_1 = -x^2 + sin(x) + 2 + rnorm(length(x),sd = 2)
y_2 = -x^2 + cos(x) - 4 + rnorm(length(x), sd = 2)
ylim = c(min(y_1, y_2), max(y_1, y_2))

y_final = c(y_2[1:200], y_1[201:400], y_2[401:length(x)])

# Time Series Plot
png("Plots/ex1_time_series.png", width = 1000, height = 600)
plot(x+3, y_final, xlab = 'Time (minutes)', ylab = "y", main = "Example Time Series",
     cex.lab=2, cex.axis=2, cex.main=2)
dev.off()
# HMM inference before FDA ----------------------------------------------------
png("Plots/ex1_hist.png", width = 700, height = 600)
hist(y_final, breaks = sqrt(length(x)), probability = T,
     xlab = "y", ylab = "", main = "Continuous Response",
     cex.lab=2, cex.axis=2, cex.main=2)

# Scaling things to plot onto the same histogram
d1 = density(y_final[201:400], adjust = 3)
scale1 = length(x) / length(201:400)
d1$y = d1$y / scale1

d2 = density(y_final[c(1:200, 401:length(x))], adjust = 3)
scale2 = length(x) / length(c(1:200, 401:length(x)))
d2$y = d2$y / scale2

lines(d1, col = 'blue', lwd = 5)
lines(d2, col = 'red', lwd=5)

dev.off()
# HMM inference with FDA ----------------------------------------------------
png("Plots/ex1_fda.png", width = 700, height = 600)
plot(x+3, y_final, xlab = 'Time (minutes)', ylab = "y", main = "Functional Response",
     cex.lab=2, cex.axis=2, cex.main=2)

lines(x+3, -x^2 + sin(x) + 2, col = 'blue', lwd = 5)
lines(x+3, -x^2 + cos(x) - 4, col = 'red', lwd = 5)

dev.off()