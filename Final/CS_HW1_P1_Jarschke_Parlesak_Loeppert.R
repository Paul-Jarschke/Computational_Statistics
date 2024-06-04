# --------------------------------------------------------------
# Computational Statistics
# Homework 1 - Problem 1
# Names: Paul Jarschke, Jan Parlesak, Leon Löppert
# --------------------------------------------------------------


# Problem 1.1: Conjugate Gradient Algorithm ----

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
  return(list(conv = conv, x = c(x)))
}


# Problem 1.2: Preconditioned Conjugate Gradient Algorithm ----

pcg <- function(A, b, x, M, tol = 1e-14) {
  # Initializations: ----
  r <- b - A %*% x  # initial residual vector
  z <- solve(M, r)  # apply preconditioner M to the residual vector
  p <- z  # direction vector
  j <- 1  # iteration counter
  
  conv <- c()  # convergence criterion tracker
  err <- c()  # error tracker
  
  conv[1] <- norm(r, "2")
  
  if (all(b == 0)) {
    # if b is a zero vector, track error
    err[1] <- norm(x, "I")
  }
  
  # Iterate until residual norm is sufficiently small
  # or reached max number of iterations
  cat("Starting Preconditioned Conguate Gradient algorithm...\n")
  
  while ((norm(r, "2") / norm(b, "2") > tol) && j < 500) {
    alpha <- crossprod(z, r) / (t(p) %*% A %*% p)  # step size
    x <- x + c(alpha) * p  # Update solution vector x
    rnew <- r - c(alpha) * (A %*% p)  # new residual vector
    znew <- solve(M, rnew)  # update z
    beta <-
      crossprod(znew, rnew) / crossprod(z, r)  # calculate beta coefficient
    p <- znew + c(beta) * p  # update direction vector
    r <- rnew  # update residual vector
    z <- znew  # update z
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
  lines(1:j, conv, lty = 2, col = "blue")
  
  plot(
    1:j,
    conv,
    type = "o",
    col = "red",
    pch = 16,
    cex = 0.5,
    main = "Convergence Plot",
    xlab = "Iteration",
    ylab = "Log Convergence Criterion"
  )
  lines(1:j, conv, lty = 2, col = "red")
  
  # Return the convergence history and the final solution vector
  return(list(conv = conv, x = c(x)))
}

