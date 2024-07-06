# --------------------------------------------------------------
# Computational Statistics
# Homework 2
# Names: Paul Jarschke, Jan Parlesak, Leon Löppert
# --------------------------------------------------------------

# ---------------------------------------------------------------
# Problem 2 - Statistical modelling using glmmTMB and R-INLA ----
# ---------------------------------------------------------------

# Initialization ----
set.seed(247)
n_obs <- 50
n_clus <- 5
n_total <- n_obs * n_clus
theta <- c(1, 2)
sigma_clus <- 1.5
lambda <- 5
cluster <- rep(1:n_clus, each=n_obs) # Cluster indicator

# Simulate data ----
x <- rnorm(n_total, 0, 1) # covariate
u <- rnorm(n_clus, 0, sigma_clus) # random effect per cluster
linear_predictor <- theta[1] + theta[2] * x + u[cluster]
lambda <- exp(linear_predictor) # lambda is the mean (prediction) of each observation
y <- rpois(n_total, lambda = lambda)
data_glmm <- data.frame(y = y, x = x, cluster = as.factor(cluster))

# 2.1) Use TMB to fit GLMM ----
# install.packages("glmmTMB)
library(glmmTMB)


# Random Intercept Model
glmm_model <-
  glmmTMB(y ~ x + (1 | cluster), family = poisson, data = data_glmm)
summary(glmm_model)

fixef(glmm_model) # Fixed Effects (FE)
ranef(glmm_model) # Random Effects (RE)
coef(glmm_model) # Total Effects (FE + RE)

# Random Intercept + Slope
glmm_model2 <-
  glmmTMB(y ~ x + (x | cluster), family = poisson, data = data_glmm)
summary(glmm_model2)

fixef(glmm_model2) # Fixed Effects (FE)
ranef(glmm_model2) # Random Effects (RE)
coef(glmm_model2) # Total Effects (FE + RE)

# 2.2) Use R-INLA to fit GLMM ----
library(INLA)

# Random Intercept Model
formula <- y ~ x + f(cluster, model = "iid")
inla_model <-
  inla(
    formula,
    family = "poisson",
    data = data_glmm
  )
summary(inla_model)

inla_model$summary.fixed # Fixed Effect distribution
inla_model$summary.random # Random Effect distribution

# Random Intercept + Slope (uncorrelated) 
cluster2 <- cluster + 5
formula2 <- y ~ x + f(cluster, model = "iid") + f(cluster2, x, model = "iid")
inla_model2 <-
  inla(
    formula2,
    family = "poisson",
    data = data_glmm
  )
summary(inla_model2)


# 2.3) Compare Results (for Random Intercept Model only) ----
# Extract fixed effect estimates
glmm_fixed <- fixef(glmm_model)
inla_fixed <- inla_model$summary.fixed

# Extract random effect estimates
glmm_random <- ranef(glmm_model)
inla_random <- inla_model$summary.random$cluster

# Print the estimates
print(glmm_fixed)
print(inla_fixed)
print(glmm_random)
print(inla_random)

# 2.4) Comparative Plot ----
