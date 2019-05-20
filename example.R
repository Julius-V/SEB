library(ggplot2)
library(gridExtra)

set.seed(3)
L <- 3
J <- 2
n <- 30
Y <- runif(n)
X1 <- runif(n)
df <- data.frame(Y = Y, X1 = X1)
p1 <- ggplot(data = df, aes(x = X1, y = Y)) + 
  geom_point() + theme_bw() + theme(panel.grid = element_blank()) +
  labs(x = expression(X[1])) + lims(x = c(0, 1), y = c(0, 1))

p2 <- p1 +
  geom_hline(yintercept = sort(Y)[11], col = 'red') +
  geom_hline(yintercept = sort(Y)[21], col = 'red')

p3 <- p2 +
  geom_segment(aes(x = sort(X1[Y < sort(Y)[11]])[6],
                   xend = sort(X1[Y < sort(Y)[11]])[6],
                   y = -Inf, yend = sort(Y)[11]), col = 'red') +
  geom_segment(aes(x = sort(X1[Y >= sort(Y)[11] & Y < sort(Y)[21]])[6],
                   xend = sort(X1[Y >= sort(Y)[11] & Y < sort(Y)[21]])[6],
                   y = sort(Y)[11], yend = sort(Y)[21]), col = 'red') +
  geom_segment(aes(x = sort(X1[Y >= sort(Y)[21]])[6],
                   xend = sort(X1[Y >= sort(Y)[21]])[6],
                   y = sort(Y)[21]), yend = Inf, col = 'red')
grid.arrange(p1, p2, p3, nrow = 1)