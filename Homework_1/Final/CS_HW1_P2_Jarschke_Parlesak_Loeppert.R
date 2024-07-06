# --------------------------------------------------------------
# Computational Statistics
# Homework 1 - Problem 2
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


# Problem 2.2: Bayesian Conjugate Gradient Method ----
bayescg <- function(A, b, x, Sig, max_it = NULL, tol = 1e-6, delay = NULL, 
                    reorth = TRUE, NormA = NULL, xTrue = NULL, 
                    SqrtSigTranspose = NULL) {
  
  ## Variable definitions: ----
  
  # Size of the system
  N <- length(x)
  
  # Default Maximum Iterations
  if (is.null(max_it)) {
    max_it <- N
  }
  
  # Residual and first search direction
  r <- matrix(0, nrow = N, ncol = max_it + 1)
  r[, 1] <- b - A(x)
  S <- r
  
  # Inner products
  rIP <- numeric(max_it + 1)
  rIP[1] <- sum(r[, 1] * r[, 1])
  sIP <- numeric(max_it)
  
  # Array holding matrix-vector products
  SigAs_hist <- matrix(0, nrow = N, ncol = max_it)
  
  # Convergence information
  # If xTrue is supplied, more information is computed
  rNorm <- sqrt(rIP[1])
  Res2 <- numeric(max_it + 1)
  if (is.null(NormA) || is.null(xTrue)) {
    bNorm <- norm(b, type = "2")
    Res <- rNorm / bNorm
    Res2[1] <- Res
  }
  
  if (!is.null(xTrue)) {
    xNorm <- norm(xTrue, type = "2")
    err_hist <- numeric(max_it + 1)
    err_hist[1] <- sum((x - xTrue) * A(x - xTrue))
    if (!is.null(NormA)) {
      xNormANorm <- norm(xTrue, type = "2") * NormA
      Res <- rNorm / xNormANorm
      Res2[1] <- Res
    }
    Res3 <- Res2
    tr_hist <- numeric(max_it + 1)
    tr_hist[1] <- sum(diag(A(Sig(diag(N)))))
  }
  
  i <- 0
  
  
  ## Iterating Through Bayesian Conjugate Gradient: ----
  
  while (i < max_it && (is.null(tol) || Res > tol)) {
    # print(paste("Iteration:", i + 1))
    # Compute Matrix Vector Products
    As <- A(S[, i + 1])
    if (!is.null(SqrtSigTranspose)) {
      SigAs_hist[, i + 1] <- SqrtSigTranspose(As)
      SigAs <- Sig(SigAs_hist[, i + 1])
    } else {
      SigAs_hist[, i + 1] <- Sig(As)
      SigAs <- SigAs_hist[, i + 1]
    }
    ASigAs <- A(SigAs)
    
    # Search Direction Inner Product
    sIP[i + 1] <- abs(sum(S[, i + 1] * ASigAs))
    
    # Calculate next x
    alpha <- rIP[i + 1] / sIP[i + 1]
    x <- x + alpha * SigAs
    
    # Calculate New Residual
    r[, i + 2] <- r[, i + 1] - alpha * ASigAs
    
    if (reorth) {
      # Reorthogonalize Residual
      r_ip_inv <- 1 / rIP[1:(i + 1)]
      ortho_term <- r[, 1:(i + 1)] %*% (t(r[, 1:(i + 1)]) %*% r[, i + 2])
      diag_r_ip_inv <- diag(r_ip_inv, nrow = length(r_ip_inv), ncol = length(r_ip_inv))
      if (ncol(ortho_term) == nrow(diag_r_ip_inv)) {
        r[, i + 2] <- r[, i + 2] - ortho_term %*% diag_r_ip_inv
      } # else {
        # if dimension checking is not successful, give warning
        # warning("Dimensions of ortho_term and diag(r_ip_inv) do not match")
      # }
    }
    
    # Compute Residual Norms
    rIP[i + 2] <- sum(r[, i + 2] * r[, i + 2])
    rNorm <- sqrt(rIP[i + 2])
    if (!is.null(xTrue)) {
      err_hist[i + 2] <- sum((x - xTrue) * A(x - xTrue))
      tr_hist[i + 2] <- tr_hist[i + 1] - sum(diag(A(outer(SigAs, SigAs)))) / sIP[i + 1]
      
      rTrueNorm <- norm(b - A(x), type = "2")
      if (!is.null(NormA)) {
        Res <- rNorm / xNormANorm
        Res3[i + 2] <- rTrueNorm / xNormANorm
      } else {
        Res3[i + 2] <- rTrueNorm / bNorm
      }
    }
    if (is.null(NormA)) {
      Res <- rNorm / bNorm
    } else if (is.null(xTrue)) {
      Res <- rNorm / NormA / norm(x, type = "2")
    }
    Res2[i + 2] <- Res
    
    # Calculate next search direction
    beta <- rIP[i + 2] / rIP[i + 1]
    S[, i + 2] <- r[, i + 2] + beta * S[, i + 1]
    
    i <- i + 1
  }
  
  ## Return results: ----
  
  info <- list(res = Res2[1:(i + 1)], search_dir = (sIP[1:i] ^ (-1/2)) * S[, 1:i])
  
  if (!is.null(xTrue)) {
    info$actual_res <- Res3[1:(i + 1)]
    info$err <- err_hist[1:(i + 1)]
    info$trace <- tr_hist[1:(i + 1)]
  }
  
  if (!is.null(delay)) {
    delay <- min(delay, i)
    post_scale <- sum(rIP[(i - delay + 1):i]^2 / sIP[(i - delay + 1):i])
    info$scale <- post_scale
  }
  
  return(list(x = x, SigAs_hist = (sIP[1:i] ^ (-1/2)) * SigAs_hist[, 1:i], info = info))
}