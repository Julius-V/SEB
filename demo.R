rm(list = ls())
library(Rcpp)
sourceCpp("source.cpp")

k <- 3 # Number of variables
C <- 3 # Number of cuts long each axis
n <- C^k * 5 # Five observations per cell
M <- matrix(rnorm(n * k, sd = 10), ncol = k)

demo1 <- SEB(M, times = rep(C, k)) # With intervals
demo2 <- SEB(M, times = rep(C, k), intervals = FALSE) # Without intervals

# Some gains!
microbenchmark::microbenchmark(SEB(M, times = rep(C, k)), SEB(M, times = rep(C, k), intervals = FALSE), times = 500)
table(demo1$labels)

# Plotting
library(RColorBrewer)
aux1 <- brewer.pal.info[brewer.pal.info$category == 'qual',]
cols <- unlist(mapply(brewer.pal, aux1$maxcolors, rownames(aux1)))
library(car)
scatter3d(M[, 1], M[, 2], M[, 3], point.col = cols[demo1$labels], surface = FALSE,
          xlab = "Y", ylab = "X1", zlab = "X2", axis.scales = FALSE, grid = FALSE, grid.lines = FALSE)

