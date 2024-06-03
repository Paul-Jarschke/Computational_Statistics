con_grad <- function(A,b,x){
  r <- b - A %*% x
  p <- r
  j <- 1
  
  conv = list()
  err = list()
  conv[1] <- norm(r, type = '2')
  
  if(all(b == 0)){
    err[1] <- norm(x, type = 'F')
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
      err[j] <- norm(x, type = 'F')
    }
  }
  print(x)
  print(j)
  print(conv[j])
  
  par(mfrow = c(1, 2))
  plot(c(1:j), conv, type = 'l', xlab = 'iteration')
  plot(c(1:j), conv, log = 'y', type ='l', xlab = 'iteration')
  
  
  
  
}


A <- array(c(4,1,1,3), dim= c(2,2))
A
b <- array(c(1,2), dim= c(2,1))
b
x_zero <- array(c(2,1), dim = c(2,1))
x_zero

con_grad(A,b,x_zero)

pre_con_grad <- function(A,b,x, M,tol){
  r <- b - A %*% x
  z <- M/r
  p <- z
  j <- 1
  
  conv = list()
  err = list()
  conv[1] <- norm(r, type = '2')
  
  if(all(b == 0)){
    err[1] <- norm(x, type = 'F')
  }
  
  while((norm(r, type = '2')/norm(b, type = '2') > tol) & j < 500){
    alpha <- (t(r) %*% r)/(t(p) %*% A %*% p)
    x <- x + as.vector(alpha) * p
    rnew <- r - as.vector(alpha) * A %*% p
    znew <- M/rnew
    beta <- (t(znew) %*% rnew)/(t(z) %*% r)
    p <- znew + as.vector(beta) * p
    r <- rnew
    z <- znew
    j <- j+1
    conv[j] <- norm(r, '2')/norm(b, '2')
    
    if(all(b == 0)){
      err[j] <- norm(x, type = 'F')
    }
  }
  print(x)
  print(j)
  print(conv[j])
  
  par(mfrow = c(1, 2))
  plot(c(1:j), conv, type = 'l', xlab = 'iteration')
  plot(c(1:j), conv, log = 'y', type ='l', xlab = 'iteration')
  
}


