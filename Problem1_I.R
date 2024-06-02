cg <- function(A, b, x) {
  # Calculate the initial residual vector (r)
  r <- b - A %*% x
  # Initialize direction vector (p) and iteration counter (j)
  p <- r
  j <- 0
  
  # Initialize convergence tracker with norm of initial residual
  conv <- numeric()
  conv[1] <- norm(r, type = "2")
  
  # If b is a zero vector, initialize error tracker
  if (all(b == 0)) {
    err <- numeric()
    err[1] <- norm(x, type = "I")
  }
  
  # Iterate until residual norm is sufficiently small or reached max number of iterations
  while ((norm(r, type = "2") / norm(b, type = "2") > 1e-14) && (j < 500)) {
    # Calculate the step size alpha
    alpha <- as.numeric((t(r) %*% r) / (t(p) %*% A %*% p))
    # Update the solution vector x
    x <- x + alpha * p
    # Calculate the new residual vector rnew
    rnew <- r - alpha * A %*% p
    # Calculate the coefficient beta
    beta <- as.numeric((t(rnew) %*% rnew) / (t(r) %*% r))
    # Update the direction vector p
    p <- rnew + beta * p
    # Update the residual vector r
    r <- rnew
    # Increment the iteration counter
    j <- j + 1
    # Track the convergence progress
    conv[j] <- norm(r, type = "2") / norm(b, type = "2")
    
    # If b is a zero vector, track the error of the solution x
    if (all(b == 0)) {
      err[j] <- norm(x, type = "I")
    }
    
    # Print the iteration number and the current convergence value
    cat(sprintf('%5.0f %20.16e\n', j, conv[j]))
  }
  
  # Set up the plotting area for two plots
  par(mfrow = c(2, 1))
  # Plot the convergence over iterations
  plot(1:j, conv, type = "o", col = "blue", pch = 16, lty = 2, xlab = "Iteration", ylab = "Convergence")
  
  # Filter out non-positive values for the logarithmic plot
  positive_indices <- which(conv > 0)
  # Plot the convergence on a logarithmic scale for positive values
  plot(positive_indices, conv[positive_indices], type = "o", log = "y", pch = 16, lty = 2, col = "red", 
       xlab = "Iteration", ylab = "Log Convergence")
  
  # Return the convergence history and the final solution vector
  return(list(conv = conv, x = x))
}





# Example usage:
set.seed(99) # For reproducibility

# Generate a random matrix
n <- 10
A <- matrix(rnorm(n * n), n, n)

# Ensure A is symmetric positive-definite
A <- t(A) %*% A + n * diag(n)

# Generate a random vector b
b <- rnorm(n)

# Initial guess for x
x0 <- rep(0, n)

result <- cg(A, b, x0)

# Print final results
cat('Final solution x:\n')
print(result$x[length(result$x)])

cat('Final convergence metric:\n')
print(result$conv[length(result$conv)])

