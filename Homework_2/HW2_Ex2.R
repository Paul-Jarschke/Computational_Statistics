# --------------------------------------------------------------
# Computational Statistics
# Homework 2
# Names: Paul Jarschke, Jan Parlesak, Leon Löppert
# --------------------------------------------------------------

# ---------------------------------------------------------------
# Problem 2 - Statistical modelling using glmmTMB and R-INLA ----
# ---------------------------------------------------------------

# Initialization ----
n_obs <- 50
n_clus <- 5
theta <- c(1, 2)
x <- rnorm(n_obs, mean = theta[1], sd = theta[2])
sigma_clus <- 1.5
cluster <- rnorm(n_clus, mean = 0, sd = sigma_clus)
lambda <- 5
y <- rpois(n_obs, lambda = lambda)


# 2.1) Use TMB to fit GLMM ----
# install.packages("glmmTMB)
library(glmmTMB)

# 2.2) Use R-INLA to fit GLMM ----
library(INLA)

# 2.3) Compare Results ----

# 2.4) Comparative Plot ----
