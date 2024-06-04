a_lanczos <- function(A, v0, max_it = NULL, tol = 1e-10, reorth = TRUE) {
  N <- length(v0)
  
  if (is.null(max_it)) {
    max_it <- N
  }
  
  V <- matrix(0, nrow = N, ncol = max_it + 1)
  
  Av <- A(v0)
  beta <- sqrt(sum(v0 * Av))
  V[, 2] <- v0 / beta  # This is index 2 intentionally, 1 not returned
  Av <- Av / beta
  
  i <- 1
  
  while (i < max_it && beta > tol) {
    # Lanczos Recursions
    w <- Av - beta * V[, i]
    alpha <- sum(w * Av)
    w <- w - alpha * V[, i + 1]
    
    # Reorthogonalization
    if (reorth) {
      w <- w - V[, 2:(i + 1)] %*% (t(V[, 2:(i + 1)]) %*% A(w))
      w <- w - V[, 2:(i + 1)] %*% (t(V[, 2:(i + 1)]) %*% A(w))
    }
    
    # Norm of New Vector
    Av <- A(w)
    beta <- sqrt(sum(w * Av))
    
    # Store New Vector
    if (beta > tol) {
      i <- i + 1
      Av <- Av / beta
      V[, i + 1] <- w / beta
    }
  }
  
  return(V[, 2:(i + 1)])
}


library(Matrix)

bayescg <- function(A, b, x, Sig, max_it = NULL, tol = 1e-6, delay = NULL, reorth = TRUE, NormA = NULL, xTrue = NULL, SqrtSigTranspose = NULL) {
  N <- length(x)
  data_type <- "double"
  
  if (is.null(max_it)) {
    max_it <- N
  }
  
  r <- matrix(0, nrow = N, ncol = max_it + 1)
  r[, 1] <- b - A(x)
  S <- r
  
  rIP <- numeric(max_it + 1)
  rIP[1] <- sum(r[, 1] * r[, 1])
  sIP <- numeric(max_it)
  
  SigAs_hist <- matrix(0, nrow = N, ncol = max_it)
  
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
  
  while (i < max_it && (is.null(tol) || Res > tol)) {
    As <- A(S[, i + 1])
    if (!is.null(SqrtSigTranspose)) {
      SigAs_hist[, i + 1] <- SqrtSigTranspose(As)
      SigAs <- Sig(SigAs_hist[, i + 1])
    } else {
      SigAs_hist[, i + 1] <- Sig(As)
      SigAs <- SigAs_hist[, i + 1]
    }
    ASigAs <- A(SigAs)
    
    sIP[i + 1] <- abs(sum(S[, i + 1] * ASigAs))
    
    alpha <- rIP[i + 1] / sIP[i + 1]
    x <- x + alpha * SigAs
    
    r[, i + 2] <- r[, i + 1] - alpha * ASigAs
    
    if (reorth) {
      r[, i + 2] <- r[, i + 2] - (rIP[1:(i + 1)]^(-1)) * r[, 1:(i + 1)] %*% (t(r[, 1:(i + 1)]) %*% r[, i + 2])
      r[, i + 2] <- r[, i + 2] - (rIP[1:(i + 1)]^(-1)) * r[, 1:(i + 1)] %*% (t(r[, 1:(i + 1)]) %*% r[, i + 2])
    }
    
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
    
    beta <- rIP[i + 2] / rIP[i + 1]
    S[, i + 2] <- r[, i + 2] + beta * S[, i + 1]
    
    i <- i + 1
  }
  
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

bayescg_random_search <- function(A, b, x, Sig, max_it = NULL, tol = 1e-6, reorth = TRUE, NormA = NULL, xTrue = NULL, SqrtSigTranspose = NULL, seed = NULL) {
  data_type <- "double"
  N <- length(x)
  
  if (is.null(max_it)) {
    max_it <- N
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  r <- b - A(x)
  inv_check <- runif(N)
  if (is.null(SqrtSigTranspose)) {
    if (norm(inv_check - A(Sig(inv_check)), type = "2") < 1e-12) {
      ASigA_func <- A
    } else {
      ASigA_func <- function(w) A(Sig(A(w)))
    }
  } else {
    if (norm(inv_check - A(Sig(SqrtSigTranspose(inv_check))), type = "2") < 1e-12) {
      ASigA_func <- A
    } else {
      ASigA_func <- function(w) A(Sig(SqrtSigTranspose(A(w))))
    }
  }
  S <- a_lanczos(ASigA_func, rnorm(N), N, 1e-15, reorth)
  max_it <- min(max_it, ncol(S))
  
  rIP <- sum(r * r)
  SigAs_hist <- matrix(0, nrow = N, ncol = max_it)
  rNorm <- sqrt(rIP)
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
  
  for (i in 1:max_it) {
    As <- A(S[, i])
    if (!is.null(SqrtSigTranspose)) {
      SigAs_hist[, i] <- SqrtSigTranspose(As)
      SigAs <- Sig(SigAs_hist[, i])
    } else {
      SigAs_hist[, i] <- Sig(As)
      SigAs <- SigAs_hist[, i]
    }
    ASigAs <- A(SigAs)
    
    sIP <- abs(sum(S[, i] * ASigAs))
    
    alpha <- rIP / sIP
    x <- x + alpha * SigAs
    
    r <- r - alpha * ASigAs
    
    rIP <- sum(r * r)
    rNorm <- sqrt(rIP)
    if (!is.null(xTrue)) {
      err_hist[i + 1] <- sum((x - xTrue) * A(x - xTrue))
      tr_hist[i + 1] <- tr_hist[i] - sum(diag(A(outer(SigAs, SigAs)))) / sIP
      
      rTrueNorm <- norm(b - A(x), type = "2")
      if (!is.null(NormA)) {
        Res <- rNorm / xNormANorm
        Res3[i + 1] <- rTrueNorm / xNormANorm
      } else {
        Res3[i + 1] <- rTrueNorm / bNorm
      }
    }
    if (is.null(NormA)) {
      Res <- rNorm / bNorm
    } else if (is.null(xTrue)) {
      Res <- rNorm / NormA / norm(x, type = "2")
    }
    Res2[i + 1] <- Res
    
    beta <- rIP / rIP
  }
  
  info <- list(res = Res2[1:(i + 1)], search_dir = (sIP^(1/2)) * S[, 1:i])
  
  if (!is.null(xTrue)) {
    info$actual_res <- Res3[1:(i + 1)]
    info$err <- err_hist[1:(i + 1)]
    info$trace <- tr_hist[1:(i + 1)]
  }
  
  return(list(x = x, SigAs_hist = (sIP^(1/2)) * SigAs_hist[, 1:i], info = info))
}