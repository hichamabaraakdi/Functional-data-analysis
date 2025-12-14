# Functional-data-analysis
mini projet
library(fda)

# Fonction de Runge
f <- function(x) {
  1 / (1 + 25 * x^2)
}

# Domaine d'étude
x_fine <- seq(-1, 1, length.out = 500)
y_fine <- f(x_fine)

# Nombre de points d'interpolation (AJOUT DE 20)
points_list <- c(3, 6, 10, 20)

# Couleurs (autant que de cas)
cols <- c("blue", "green", "purple", "orange")

# Plot de la fonction exacte
plot(x_fine, y_fine, type="l", lwd=3, col="black",
     xlab="x", ylab="f(x)",
     main="Approximation de la fonction de Runge par B-splines")
grid()

for (i in seq_along(points_list)) {
  
  n <- points_list[i]
  
  # Points d'interpolation
  x_nodes <- seq(-1, 1, length.out = n)
  y_nodes <- f(x_nodes)
  
  # Nombre de bases = nombre de points
  nbasis <- n
  
  # Ordre des splines (linéaire)
  norder <- 2
  
  # Création de la base spline
  basis_bspline <- create.bspline.basis(rangeval = c(-1, 1),
                                        nbasis = nbasis,
                                        norder = norder)
  
  # Approximation spline
  fd_obj <- smooth.basis(argvals = x_nodes, y = y_nodes,
                         fdParobj = basis_bspline)$fd
  
  # Évaluation sur x_fine
  y_spline <- eval.fd(x_fine, fd_obj)
  
  # Tracé
  lines(x_fine, y_spline, col = cols[i], lwd = 2)
  points(x_nodes, y_nodes, col = cols[i], pch = 19)
}

legend("topright",
       legend = c("f(x) exacte",
                  paste("Approximation avec", points_list, "points")),
       col = c("black", cols),
       lty = 1,
       lwd = c(3, rep(2, length(points_list))),
       pch = c(NA, rep(19, length(points_list))))
