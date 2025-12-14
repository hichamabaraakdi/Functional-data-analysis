library(fda)
library(splines)

set.seed(123)
BM <- cumsum(rnorm(10000))
idx <- seq(1, 10000, length.out = 100)
y <- BM[idx]
n <- length(y)
t <- seq(0, 1, length.out = n)

p <- 25
degree <- 3
Phi <- bs(t, degree = degree, df = p, intercept = TRUE)

# Création de la matrice R
D2 <- diff(Phi, differences = 2)  
R <- t(D2) %*% D2

lambda_grid <- 10^seq(-12, 12, length.out = 800)
PhiTPhi <- t(Phi) %*% Phi
PhiTy   <- t(Phi) %*% y

res_norm <- numeric(length(lambda_grid))
pen_norm <- numeric(length(lambda_grid))
C_hat    <- matrix(0, nrow = length(lambda_grid), ncol = p)

#Résolution
for (i in seq_along(lambda_grid)) {
  lambda <- lambda_grid[i]
  c_hat <- solve(PhiTPhi + lambda * R, PhiTy)
  C_hat[i, ] <- c_hat
  
  # ||Phi c - y||
  res_norm[i] <- sqrt(sum((Phi %*% c_hat - y)^2))
  
  # sqrt(c^T R c)
  pen_norm[i] <- sqrt(as.numeric(t(c_hat) %*% R %*% c_hat))
}

#Plot
plot(res_norm, pen_norm,
     log = "xy",
     type = "l",
     lwd = 2,
     xlab = expression('||Phi . c - y||'),
     ylab = expression('sqrt(c^T . R . c)'),
     main = "L-curve (B-splines cubiques)")

#Calcul de la courbure 
log_r <- log(res_norm)
log_p <- log(pen_norm)

d1_r <- diff(log_r)
d1_p <- diff(log_p)

d2_r <- diff(d1_r)
d2_p <- diff(d1_p)

curvature <- abs(d1_r[-1] * d2_p - d1_p[-1] * d2_r) /
  ( (d1_r[-1]^2 + d1_p[-1]^2)^(3/2) )

idx_opt <- which.max(curvature) + 1
lambda_opt <- lambda_grid[idx_opt]

points(res_norm[idx_opt], pen_norm[idx_opt],
       col = "red", pch = 19, cex = 1.3)

cat("Lambda optimal (L-curve) =", lambda_opt, "\n")


#Courbe
c_opt <- C_hat[idx_opt, ]

plot(t, y,
     pch = 16, cex = 0.6,
     xlab = "t", ylab = "y")

lines(t, Phi %*% c_opt,
      col = "red", lwd = 2)
