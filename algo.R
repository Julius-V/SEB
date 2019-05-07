rm(list = ls())
xOHC <- function(X, C) {
  
}
partOHC <- function(Y, X, L, C) {
  n <- nrow(X)
  thrs <- sort(Y)[floor(n / L * (1:L))]
  nums0 <- as.numeric(cut(Y, c(-Inf, thrs), labels = 1:L, include.lowest = TRUE))
  maxs <- apply(X, 1, max)
  subX <- lapply(split(X, nums0), matrix, ncol = ncol(X))
  nums <- numeric(n)
  for(i in 1:L)
    nums[nums0 == i] <- xOHC()
  subX <- lapply(subX, xOHC, C)
}

partSEB <- function(D, L = NULL, C) {
  nr <- nrow(D)
  nc <- ncol(D)
  col <- D[, 1]
  nints <- if (is.null(L)) C else L
  thrs <- sort(col)[floor(nr / nints * (1:nints))]
  nums0 <- as.numeric(cut(col, c(-Inf, thrs), labels = 1:nints, include.lowest = TRUE))
  if(nc > 1) {
    ms <- lapply(split(D[, -1, drop = FALSE], nums0), matrix, ncol = nc - 1)
    nums <- numeric(nr)
    for(i in 1:nints)
      nums[nums0 == i] <- partSEB(ms[[i]], C = C) + C^(nc - 1) * (i - 1)
    nums
  } else {
    nums0
  }
}
nc <- 3
L <- 5
C <- 3
# D <- matrix(c(runif(100000), rnorm(100000), rlnorm(100000, sdlog = 0.5)), ncol = nc)
D <- matrix(runif(50000 * nc), ncol = nc)
table(partSEB(D, L, C))
L * C^(nc - 1)

library(RColorBrewer)
n <- L * C^(nc - 1)
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
library(car)
scatter3d(D[, 1], D[, 2], D[, 3], point.col = col_vector[partSEB(D, L, C)], surface = FALSE)
# plot(D[, 2], partSEB(D, L, C))

