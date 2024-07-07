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
theta <- c(1, 0.5)
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

# Plot data ----
library(ggplot2)
data_plot <-
  ggplot(data=data_glmm, aes(x, y, color = cluster)) +
  geom_point(alpha=0.4) +
  labs(title = "Data") +
  theme(plot.title = element_text(hjust = 0.5))
data_plot

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

# Plot fitted values
data_plot +
  geom_line(aes(x, fitted(glmm_model), color = cluster), linetype = "dashed")

# # Random Intercept + Slope (ignore)
# glmm_model2 <-
#   glmmTMB(y ~ x + (x | cluster), family = poisson, data = data_glmm)
# summary(glmm_model2)
# 
# fixef(glmm_model2) # Fixed Effects (FE)
# ranef(glmm_model2) # Random Effects (RE)
# coef(glmm_model2) # Total Effects (FE + RE)

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

# # Random Intercept + Slope (uncorrelated) (ignore)
# cluster2 <- cluster + 5
# formula2 <- y ~ x + f(cluster, model = "iid") + f(cluster2, x, model = "iid")
# inla_model2 <-
#   inla(
#     formula2,
#     family = "poisson",
#     data = data_glmm,
#     control.predictor = list(compute = TRUE),
#     control.compute=list(return.marginals.predictor=TRUE)
#   )
# summary(inla_model2)

# Plot fitted values
data_plot +
  geom_line(aes(x, c(inla_model$summary.fitted.values["mean"])$mean, color = cluster), linetype = "dotted")

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

# 2.4) Comparative Plots ----
library(ggplot2)

# Compare fitted values between TMB and INLA
glmm_fitted <- data.frame(
  x = x,
  TMB = fitted(glmm_model),
  INLA = c(inla_model$summary.fitted.values["mean"])$mean,
  cluster = data_glmm$cluster
)

data_plot +
  geom_line(data = glmm_fitted, aes(x = x, y = TMB, linetype = "TMB")) +
  geom_line(data = glmm_fitted, aes(x = x, y = INLA, linetype = "INLA")) +
  scale_linetype_manual(name = "Fit", values = c("dashed", "dotted"), 
                        breaks = c("TMB", "INLA"), labels = c("TMB", "INLA")) +
  labs(title = "Fitted Models Comparison", x = "X Axis", y = "Y Axis") 
  # identical mean predictions

# Compare fixed effect estimates
fixed_effects <- data.frame(
  Model = c("glmmTMB", "glmmTMB", "R-INLA", "R-INLA"),
  Parameter = rep(c("Intercept", "Slope"), 2),
  Estimate = c(glmm_fixed$cond[1], glmm_fixed$cond[2], inla_fixed$mean[1], inla_fixed$mean[2]),
  Lower = c(ci_glmm[1, 1], ci_glmm[2, 1], inla_fixed$`0.025quant`[1], inla_fixed$`0.025quant`[2]),
  Upper = c(ci_glmm[1, 2], ci_glmm[2, 2], inla_fixed$`0.975quant`[1], inla_fixed$`0.975quant`[2])
)

# compare mode of INLA estimation not mean to fixed effects estimate from TMB?

ggplot(fixed_effects, aes(x = Parameter, y = Estimate, color = Model)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, position = position_dodge(width = 0.5)) +
  labs(title = "Fixed Effect Estimates (+ 95% CI)", y = "Estimate")

# You can not compare random effects from TMB and R-INLA graphically in the same way as above,
# since TMB does not provide estimates for uncertainty of Random Effects 
# as it does for the fixed effects. 
# R-INLA, however, does so:
random_effects_df <-
  data.frame(
    cluster = inla_random$ID,
    Mean_INLA = inla_random$mean,
    Lower = inla_random$`0.025quant`,
    Upper = inla_random$`0.975quant`,
    Mean_TMB = c(glmm_random$cond$cluster)$`(Intercept)`
  )
library(ggplot2)

# Creating the ggplot
ggplot(random_effects_df, aes(x = cluster)) +
  geom_point(aes(y = Mean_INLA, color = "INLA", shape = "INLA"), size = 1.2) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = "INLA"), width = 0.2) +
  geom_point(aes(y = Mean_TMB, color = "TMB", shape = "TMB"), size = 0.8) +
  labs(title = "Random Effect Estimates (R-INLA)",
       y = "Estimate",
       color = "Legend",
       shape = "Legend") +
  scale_color_manual(values = c("INLA" = "black", "TMB" = "red")) +
  scale_shape_manual(values = c("INLA" = 16, "TMB" = 16)) +
  theme(legend.position = "right")
  # this graph is not really pretty but you can see that the mean estimates are similar

# Plot Posterior Distribution of R-INLA Parameter Estimates
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