# --------------------------------------------------------------
# Computational Statistics
# Homework 1
# Name: Leon Löppert
# --------------------------------------------------------------

# Load dependencies ----
library(pracma)

# Problem 1.2 Preconditioned Conjugate Gradient Algorithm ----

PCG <-
  function(A, b, x, M, tol = 1e-14) {
    r <- b - A %*% x
    z <- solve(M, r)
    p <- z
    j <- 1
    
    conv <- Norm(r)
    
    if(all(b == 0)) {
      error <- Norm(r, Inf)
    }
    
    cat("Starting Preconditioned Conguate Gradient algorithm...\n")
    while (((Norm(r) / Norm(b)) > tol) && j < 500) {
      alpha <- crossprod(z, r) / (t(p) %*% A %*% p)
      x <- x + c(alpha) * p
      rnew <- r - c(alpha) * (A %*% p)
      znew <- solve(M, rnew)
      beta <- crossprod(znew, rnew) / crossprod(z, r)
      p <- znew + c(beta) * p
      r <- rnew
      z <- znew
      j <- j + 1
      
      conv[j] <- Norm(r) / Norm(b)
      
      if(all(b == 0)) {
        err[j] <- Norm(x, Inf)
      }
      cat(sprintf('It. %3.0f: %20.16e', j, conv[j]), "\n")
    }
    cat("...finished!\n")
    cat("The solution for x is:", x, "\n")
    
    par(mfrow = c(2, 1))
    
    plot(
      1:j,
      conv,
      type = 'l',
      col = 'black',
      main = 'Convergence Plot',
      xlab = 'Iteration',
      ylab = 'Conv. Criterion', cex = 0.5
    )
    lines(1:j, conv, lty = 2, col = 'black')
    
    plot(
      1:j,
      conv,
      type = 'o',
      log = 'y',
      col = 'black',
      main = 'Convergence Plot',
      xlab = 'Iteration',
      ylab = 'Conv. Criterion (log-scale)',
      cex = 0.5
    )
    lines(1:j, conv, lty = 2, col = 'black')
    
    return(list(conv = conv, x = x))
  }


# Testing algorithm ----
n <- 5
generateSPDMatrix <- function(n) {
  A <- matrix(rnorm(n * n), n, n)
  A <- A %*% t(A)  # Make it symmetric positive definite
  A <- A + n * diag(n)  # Ensure positive definiteness
  return(A)
}

A <- generateSPDMatrix(n)
b <- rnorm(n)
x <- 1:n
b <- A %*% x

x0 <- rep(0, n)
M <- diag(n)  # Preconditioner (here: identity matrix)

PCG_results <- PCG(A, b, x0, M)
PCG_results

R_solver_results <- solve(A, b)
R_solver_results
