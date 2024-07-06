# --------------------------------------------------------------
# Computational Statistics
# Homework 2
# Names: Paul Jarschke, Jan Parlesak, Leon Löppert
# --------------------------------------------------------------

# -----------------------------
# Problem 1 - EM Algorithm ----
# -----------------------------

# 1.1) ----

# Set seed for reproducibility
set.seed(1234)
n <- 2000

# Generate random uniform variable (covariate x)
x <- runif(n, -5, 5)

# Initialize parameters
beta <- c(0, 1.5)
sigma <- 1.0
a <- 20
b <- a
pi <- 0.45
theta <- c(beta, sigma, pi)
c <- 1/(2*a)  # uniform distribution is symmetric around zero


# 1.2) ----

# Define function to generate y
generate_y <- function(x, beta, sigma, a, b, pi) {
  n <- length(x)
  y <- numeric(n)
  for (i in 1:n) {
    if (runif(1) < pi) {
      y[i] <- rnorm(1, mean = beta[1] + beta[2] * x[i], sd = sigma)
    } else {
      y[i] <- runif(1, -a, b)
    }
  }
  return(y)
}


# Generate response y
y <- generate_y(x, beta, sigma, a, b, pi)

# Store data in data frame
data <- data.frame(x = x, y = y)
head(data)

# Fit OLS regression model
ols_model <- lm(y ~ x, data = data)
coef(ols_model)

# Plot data and regression line
plot(
  data$x,
  data$y,
  main = "OLS Regression",
  xlab = "x",
  ylab = "y",
  pch = 19,
  col = "blue"
)
abline(ols_model, col = "red", lwd = 2)
legend("topright", legend = c("Data", "Fit"), col = c("blue", "red"), 
       pch = c(19, NA), lty = c(NA, 1), lwd = c(NA, 2))



############ (3) ############

em_algorithm <- function(data, theta_init, a, max_iter = 100, tol = 1e-6) {
  n <- nrow(data)
  x <- data$x
  y <- data$y
  
  beta <- theta_init[[1]]
  sigma <- theta_init[[2]]
  pi <- theta_init[[3]]

  for (iter in 1:max_iter) {
    # E-step: Calculate posterior probabilities
    p1 <- pi * dnorm(y, mean = beta[1] + beta[2] * x, sd = sigma)
    p2 <- (1 - pi) * dunif(y, -a, b)
    Ez <- p1 / (p1 + p2)

    # M-step: Update parameters
    lm_fit <- lm(y ~ x, weights = Ez)
    beta_new <- coef(lm_fit)
    sigma_new <- sqrt(sum(Ez * (y - predict(lm_fit))^2) / sum(Ez))
    pi_new <- mean(Ez)

    # Check for convergence
    if (max(abs(beta_new - beta), abs(sigma_new - sigma), abs(pi_new - pi)) < tol) {
      break # break when all parameters are converged
    }
    # Update variables
    beta <- beta_new
    sigma <- sigma_new
    pi <- pi_new
  }
  # Create result list
  list(beta = beta, sigma = sigma, pi = pi, iter = iter)
}

# Run EM algorithm
beta_init <- c(0, 1)
sigma_init <- 4
pi_init <- 0.5
theta_init <- list(beta_init, sigma_init, pi_init)

em_result <- em_algorithm(data, theta_init, a)
em_result
