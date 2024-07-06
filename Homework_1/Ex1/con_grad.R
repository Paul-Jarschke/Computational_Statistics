# --------------------------------------------------------------
# Computational Statistics
# Homework 1
# Names: Paul Jarschke, Jan Parlesak, Leon Löppert
# --------------------------------------------------------------

# Problem 1.1 Conjugate Gradient Algorithm ----

cg <- function(A, b, x) {
  # Initializations:
  r <- b - A %*% x  # initial residual vector
  p <- r  # direction vector
  j <- 0  # iteration counter
  
  conv <- c()  # convergence criterion tracker
  err <- c()  # error tracker
  
  conv[1] <- norm(r, "2")
  
  if (all(b == 0)) {
    # if b is a zero vector, track error
    err[1] <- norm(x, "I")
  }
  
  # Iterate until residual norm is sufficiently small
  # or reached max number of iterations
  cat("Starting Conguate Gradient algorithm...\n")
  while ((norm(r, "2") / norm(b, "2") > 1e-14) && j < 500) {
    alpha <- (crossprod(r)) / (t(p) %*% A %*% p)  # step size
    x <- x + c(alpha) * p  # Update solution vector x
    rnew <- r - c(alpha) * (A %*% p)  # new residual vector
    beta <-
      crossprod(rnew) / crossprod(r)  # calculate beta coefficient
    p <- rnew + c(beta) * p  # update direction vector
    r <- rnew  # update residual vector
    j <- j + 1  # update iteration counter
    conv[j] <- norm(r, "2") / norm(b, "2")  # track convergence
    
    if (all(b == 0)) {
      # if b is a zero vector, track error
      err[j] <- norm(x, "I")
    }
    # Iteration counter:
    cat(sprintf('It. %3.0f: %20.16e\n', j, conv[j]), "\n")
  }
  cat("...finished!\n")
  cat("The solution for x is:", x, "\n")
  
  # Convergence plots (Convergence criterion over iterations):
  par(mfrow = c(2, 1))
  
  plot(
    1:j,
    conv,
    type = "o",
    col = "blue",
    pch = 16,
    cex = 0.5,
    main = "Convergence Plot",
    xlab = "Iteration",
    ylab = "Convergence Criterion"
  )
  lines(1:j, conv, lty = 2, col = 'blue')
  
  plot(
    1:j,
    conv,
    log = "y",
    type = "o",
    col = "red",
    pch = 16,
    cex = 0.5,
    main = "Convergence Plot",
    xlab = "Iteration",
    ylab = "Log Convergence Criterion"
  )
  lines(1:j, conv, lty = 2, col = 'red')
  
  # Return the convergence history and the final solution vector
  return(list(conv = conv, x = x))
}

# Testing CG algorithm ----

## Easy example: ----
A <- matrix(c(4, 1, 1, 3), nrow = 2)  # Coefficient matrix A
x <- matrix(c(1, 2), nrow = 2)  # True x
b <- A %*% x  # b vector (RHS)

x0 <- c(0, 0)  # initial guess for x

cg_results <- cg(A, b, x0)
R_solver_results <- solve(A, b)

x
cg_results
R_solver_results

## More complex examples: ----

# Generate a random matrix
n <- 10
A <- matrix(rnorm(n * n), n, n)

# Ensure A is symmetric positive-definite
A <- crossprod(A) + n * diag(n)

# Generate a random vector b
b <- rnorm(n)

# Initial guess for x
x0 <- rep(0, n)

cg_results <- cg(A, b, x0)
R_solver_results <- solve(A, b)

# Print final results
cat('Final solution x:\n')
print(cg_results$x)
print(R_solver_results)