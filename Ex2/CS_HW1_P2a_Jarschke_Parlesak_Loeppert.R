# --------------------------------------------------------------
# Computational Statistics
# Homework 1
# Name: Paul Jarschke, Jan Parlesak, Leon Loeppert
# --------------------------------------------------------------

# Problem 2.1: Modified Lanczos that creates A-orthonormal bases ----

a_lanczos <-
  function(A,  # Input matrix
           v0,  # Initial vector
           max_it = NULL,  # Maximum number of iterations
           tol = 1e-10,  # Tolerance value
           reorth = TRUE) {  # Flag for reorthogonalization
    
    # Preliminary calculations and initializations: ----
    N <- length(v0)
    
    if(!is.null(max_it)) max_it <- N
    
    # Initialize the basis matrix V with zeros
    V <- matrix(0, nrow = N, ncol = max_it + 1)
    
    # Compute initial vector Av and its norm beta
    Av <- A(v0)
    beta <- sqrt(sum(v0 * Av))
    
    # Set the first basis vector
    V[, 2] <-
      v0 / beta  # This is index 2 intentionally, 1 not returned
    
    # Normalize the initial vector Av
    Av <- Av / beta
    
    i <- 2  # Iteration counter
    
    # Lanczos iteration loop
    while (i < (max_it + 1) && beta > tol) {
      # Lanczos recursions
      w <- Av - beta * V[, i - 1]
      alpha <- sum(w * Av)
      w <- w - alpha * V[, i]
      
      # Reorthogonalization step
      if (reorth) {
        # Subtract projections onto previous basis vectors (twice)
        w <- w - (V[, 2:i] %*% crossprod(V[, 2:i], A(w)))
        w <- w - (V[, 2:i] %*% crossprod(V[, 2:i], A(w)))
      }
      
      # Norm of New Vector Av
      Av <- A(w)
      beta <- sqrt(sum(w * Av))
      
      # Store new vector if its norm is above the tolerance
      if (beta > tol) {
        i <- i + 1
        Av <- Av / beta
        V[, i] <- w / beta
      }
    }
    
    # Return the matrix V containing the A-orthonormal basis vectors
    return(V[, 2:i])
  }