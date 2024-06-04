pcg <- function(A, b, x, M, tol) {
  # Calculate the initial residual vector r = b - A*x
  r <- b - A %*% x
  # Apply the preconditioner M to the residual vector
  z <- solve(M, r)
  p <- z
  j <- 1
  
  # Initialize the convergence tracker
  conv <- numeric()
  conv[1] <- norm(r)
  
  # If b is a zero vector, initialize the error tracker
  if (all(b == 0)) {
    err <- numeric()
    err[1] <- norm(x, type = "I")
  }
  
  # Iterate until the residual norm is sufficiently small or the maximum number of iterations is reached
  while ((norm(r) / norm(b) > tol) && (j < 500)) {
    # Calculate the step size alpha
    alpha <- t(z) %*% r / (t(p) %*% A %*% p)
    # Update the solution vector x
    x <- x + alpha * p
    # Calculate the new residual vector rnew
    rnew <- r - alpha * A %*% p
    # Apply the preconditioner M to the new residual vector
    znew <- solve(M, rnew)
    # Calculate the coefficient beta
    beta <- t(znew) %*% rnew / (t(z) %*% r)
    # Update the direction vector p
    p <- znew + beta * p
    # Update the residual and preconditioned vectors
    r <- rnew
    z <- znew
    j <- j + 1
    # Track the convergence progress
    conv[j] <- norm(r) / norm(b)
    
    # If b is a zero vector, track the error of the solution x
    if (all(b == 0)) {
      err[j] <- norm(x, type = "I")
    }
  }
  
  # Plot the convergence over iterations with circles and dashed lines
  plot(1:j, conv, type = "o", col = "blue", pch = 16, lty = 2, xlab = "Iteration", ylab = "Convergence")
  # Plot the convergence on a logarithmic scale with circles and dashed lines
  plot(1:j, conv, type = "o", col = "blue", pch = 16, lty = 2, log = "y", xlab = "Iteration", ylab = "Log Convergence")
  
  # Return the convergence history and the final solution vector
  return(list(conv = conv, x = x))
  
  
}
  







  
