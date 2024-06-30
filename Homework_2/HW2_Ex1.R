############ (1) ############

# Set seed for reproducibility
set.seed(1234)
n <- 200    #much better results/estimation with 2000 or even 20000 (longer)
x <- runif(n, -5, 5)



############ (2) ############

# Define the function to generate y
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

# Define the parameters
beta <- c(0, 1.5)
sigma <- 1.0
a <- 20
b <- 20
pi <- 0.45

# Generate the response y using the defined function
y <- generate_y(x, beta, sigma, a, b, pi)

# Create a data frame with the generated data
data <- data.frame(x = x, y = y)

# Look at data
head(data)

# Fit OLS regression model
ols_model <- lm(y ~ x, data = data)

# Plot the data and the regression line
plot(data$x, data$y, main = "OLS Regression", xlab = "x", ylab = "y", pch = 19, col = "blue")
abline(ols_model, col = "red", lwd = 2)




############ (3) ############

em_algorithm <- function(data, beta, sigma, a, pi, max_iter = 100, tol = 1e-6) {
  n <- nrow(data)
  x <- data$x
  y <- data$y

  for (iter in 1:max_iter) {
    # E-step: Calculate the posterior probabilities
    p1 <- pi * dnorm(y, mean = beta[1] + beta[2] * x, sd = sigma)
    p2 <- (1 - pi) * dunif(y, -a, b)
    w <- p1 / (p1 + p2)

    # M-step: Update the parameters
    lm_fit <- lm(y ~ x, weights = w)
    beta_new <- coef(lm_fit)
    sigma_new <- sqrt(sum(w * (y - predict(lm_fit))^2) / sum(w))
    pi_new <- mean(w)

    # Check for convergence
    if (max(abs(beta_new - beta), abs(sigma_new - sigma), abs(pi_new - pi)) < tol) {
      break
    }

    beta <- beta_new
    sigma <- sigma_new
    pi <- pi_new
  }

  list(beta = beta, sigma = sigma, pi = pi, iter = iter)
}

# Initial values
beta_init <- c(0, 1.5)
sigma_init <- 1.0
pi_init <- 0.45

# Run the EM algorithm
em_result <- em_algorithm(data, beta_init, sigma_init, a, pi_init)
em_result




















