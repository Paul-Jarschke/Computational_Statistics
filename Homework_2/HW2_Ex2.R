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
    data = data_glmm,
    control.predictor = list(compute = TRUE),
    control.compute=list(return.marginals.predictor=TRUE)
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
    data = data_glmm,
    control.predictor = list(compute = TRUE),
    control.compute=list(return.marginals.predictor=TRUE)
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

ci_glmm <- confint(glmm_model)
print(ci_glmm)

# 2.4) Comparative Plot ----
# Plot fixed effect estimates
fixed_effects <- data.frame(
  Model = c("glmmTMB", "glmmTMB", "R-INLA", "R-INLA"),
  Parameter = rep(c("Intercept", "Slope"), 2),
  Estimate = c(glmm_fixed$cond[1], glmm_fixed$cond[2], inla_fixed$mean[1], inla_fixed$mean[2]),
  Lower = c(ci_glmm[1, 1], ci_glmm[2, 1], inla_fixed$`0.025quant`[1], inla_fixed$`0.025quant`[2]),
  Upper = c(ci_glmm[1, 2], ci_glmm[2, 2], inla_fixed$`0.975quant`[1], inla_fixed$`0.975quant`[2])
)

library(ggplot2)
ggplot(fixed_effects, aes(x = Parameter, y = Estimate, color = Model)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, position = position_dodge(width = 0.5)) +
  labs(title = "Fixed Effect Estimates (+ 95% CI)", y = "Estimate")

# Plot random effect estimates
inla_random_df <- data.frame(ID = inla_random$ID, Mean = inla_random$mean, Lower = inla_random$`0.025quant`, Upper = inla_random$`0.975quant`)
ggplot(inla_random_df, aes(x = ID, y = Mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  labs(title = "Random Effect Estimates (R-INLA)", y = "Estimate")

# Plot Posterior Distribution of R-Inla Parameter Estimates
marginals_fixed <- inla_model$marginals.fixed
marginals_random <- inla_model$marginals.random$cluster
names(marginals_random) <- paste0("random_intercept.", 1:n_clus)
list_marginals <- c(marginals_fixed, marginals_random)
list_marginals

marginals <- data.frame(do.call(rbind, list_marginals))
marginals$parameter <- rep(names(list_marginals),
                          times = sapply(list_marginals, nrow))

ggplot(marginals, aes(x = x, y = y)) + geom_line() +
  facet_wrap(~ parameter,scales = "free") +
  labs(x = "", y = "Density", title = "Posterior Distribution of Parameters") +
  theme_bw()

