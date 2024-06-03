# --------------------------------------------------------------
# Computational Statistics
# Homework 1
# Name: Leon Löppert
# --------------------------------------------------------------

# Load dependencies ----
library(pracma)

# Problem 1.1 Conjugate Gradient Algorithm ----

CG <-
  function(A, b, x) {
    r <- b - A %*% x  # Initialize residual vector
    p <- r
    j <- 0
    
    conv <- Norm(r)
    
    if(all(b == 0)) {
      error <- Norm(r, Inf)
    }
    
    cat("Starting Conguate Gradient algorithm...\n")
    while (((Norm(r) / Norm(b)) > 1e-14) && j < 500) {
      alpha <- crossprod(r) / (t(p) %*% A %*% p)
      x <- x + c(alpha) * p
      rnew <- r - c(alpha) * (A %*% p)
      beta <- crossprod(rnew) / crossprod(r)
      p <- rnew + c(beta) * p
      r <- rnew
      j <- j + 1
      
      conv[j] <- Norm(r) / Norm(b)
      
      if(all(b == 0)) {
        err[j] <- Norm(x, Inf)
      }
      cat(sprintf('It. %3.0f: %20.16e', j, conv[j]), "\n")
    }
    cat("...finished! ")
    cat("The solution for x is:", x, "\n")

    par(mfrow = c(2, 1))
    
    plot(
      1:j,
      conv,
      type = 'l',
      col = 'black',
      main = 'Convergence Plot',
      xlab = 'Iteration',
      ylab = 'Conv. criterion', cex = 0.5
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
      ylab = 'Conv. criterion (log-scale)',
      cex = 0.5
    )
    lines(1:j, conv, lty = 2, col = 'black')
    
    return(list(conv = conv, x = x))
  }


# Testing algorithm ----
A <- matrix(c(4, 1, 2, 3), nrow = 2, byrow = TRUE)
x <- c(1, 2)
b <- A %*% x

x0 <- c(0,0)

CG_results <- CG(A, b, x0)
#CG_results

R_solver_results <- solve(A, b)
R_solver_results



