con_grad <- function(A,b,x){
  r <- b - A %*% x
  p <- r
  j <- 0
  
  conv = list()
  err = list()
  conv[1] <- norm(r, type = '2')
  
  if(all(b == 0)){
    err[1] <- norm(x, type = 'I')
  }
  
  while((norm(r, type = '2')/norm(b, type = '2') > 1.e-14) & j < 500){
    alpha <- (t(r) %*% r)/(t(p) %*% A %*% p)
    x <- x + as.vector(alpha) * p
    rnew <- r - as.vector(alpha) * A %*% p
    beta <- (t(rnew) %*% rnew)/(t(r) %*% r)
    p <- rnew + as.vector(beta) * p
    r <- rnew
    j <- j+1
    conv[j] <- norm(r, '2')/norm(b, '2')
    
    if(all(b == 0)){
      err[j] <- norm(x, type = 'I')
    }
  }
  cat(sprintf('%5.0f %20.16e\n', j, conv[j]))
  
  par(mfrow = c(2, 1))
  
  plot(1:j, conv, type = "o", col = "blue", pch = 16, lty = 2, xlab = "Iteration", ylab = "Convergence")
  
  positive_indices <- which(conv > 0)
  plot(positive_indices, conv[positive_indices], type = "o", log = "y", pch = 16, lty = 2, col = "red", 
       xlab = "Iteration", ylab = "Log Convergence")
  
  # Return the convergence history and the final solution vector
  return(list(conv = conv, x = x))
  

  
}


A <- array(c(4,1,1,3), dim= c(2,2))
A
b <- array(c(1,2), dim= c(2,1))
b
x_zero <- array(c(2,1), dim = c(2,1))
x_zero

con_grad(A,b,x_zero)


pre_con_grad <- function(A,b,x,M,tol){
  r <- b - A %*% x
  z <- solve(M, r)
  p <- z
  j <- 1
  
  conv = list()
  err = list()
  conv[1] <- norm(r, type = '2') 
  
  if(all(b == 0)){
    err[1] <- norm(x, type = 'I')
  }
  
  while((norm(r, type = '2')/norm(b, type = '2') > tol) && j < 500){
    alpha <- (t(z) %*% r)/(t(p) %*% A %*% p)
    x <- x + as.vector(alpha) * p
    rnew <- r - as.vector(alpha) * A %*% p
    znew <- solve(M, rnew)
    beta <- (t(znew) %*% rnew)/(t(z) %*% r)
    p <- znew + as.vector(beta) * p
    r <- rnew
    z <- znew
    j <- j+1
    conv[j] <- norm(r, '2')/norm(b, '2')
    
    if(all(b == 0)){
      err[j] <- norm(x, type = 'I')
    }
  }
  cat(sprintf('%5.0f %20.16e\n', j, conv[j]))
  
  plot(1:j, conv, type = "o", col = "blue", pch = 16, lty = 2, xlab = "Iteration", ylab = "Convergence")
  plot(1:j, conv, type = "o", col = "blue", pch = 16, lty = 2, log = "y", xlab = "Iteration", ylab = "Log Convergence")
  
  return(list(conv = conv, x = x))
  
}

M <- diag(c(4,3)) #Jacobi preconditioner 
#https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://www.cse.psu.edu/~b58/cse456/lecture20.pdf&ved=2ahUKEwi99f28oMCGAxVo2gIHHT67AY4QFnoECBYQAQ&usg=AOvVaw1PLAHorVnnAUS7vzuSUq69
M

pre_con_grad(A,b,x_zero,M, 1.e-14)
