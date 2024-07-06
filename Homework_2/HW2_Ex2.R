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
glmm_model <- glmmTMB(y ~ x + (1 | cluster), family = poisson, data = data_glmm)
summary(glmm_model)
coef(glmm_model)
ranef(glmm_model)

# Mixed Effects Model (Random Intercept + Slope)
glmm_model2 <- glmmTMB(y ~ x + (x | cluster), family = poisson, data = data_glmm)
summary(glmm_model2)
coef(glmm_model2)
ranef(glmm_model2)

# 2.2) Use R-INLA to fit GLMM ----
library(INLA)

# 2.3) Compare Results ----

# 2.4) Comparative Plot ----
