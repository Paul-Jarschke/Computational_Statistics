# --------------------------------------------------------------
# Computational Statistics
# Homework 1 - Application and Examples
# Names: Paul Jarschke, Jan Parlesak, Leon Loeppert
# --------------------------------------------------------------

# Source Code
sapply(list.files(pattern = "^CS_HW1_P", recursive = TRUE), source)

# Set seed for reproducability
set.seed(123)

# Number of rows in linear system of equations
n <- 10

# Create a random symmetric positive-definite coefficient matrix A
createPSDmatrixA <-
  function(n) {
    A <- matrix(rnorm(n * n), n, n)
    A <- crossprod(A) + n * diag(n)
  }
A <- createPSDmatrixA(n)

# Generate a random vector b
b <- rnorm(n)

# Initial guess for x
x0 <- rep(0, n)

# Find solution using Conjugate Gradient Algorithm (Problem 1.1) ----
cg_results <- cg(A, b, x0)

# Find solution using Preconditioned Conjugate Gradient Algorithm (Problem 1.2) ----
M1 <- diag(diag(A))  # Jacobi preconditioner
M2 <- diag(n)  # Identity matrix

pcg_results1 <- pcg(A, b, x0, M1)
pcg_results2 <- pcg(A, b, x0, M2)

# Testing a_lanczos function (Problem 2.1) ----
v0 <- rnorm(n)
lanczos_vectors <- a_lanczos(A_func, v0, max_it = 10, tol = 1e-10, reorth = TRUE)
print("Lanczos vectors:")
lanczos_vectors

# Find solution using the Bayesian Conjugate Gradient Algorithm (Problem 2.2) ----
Sig <- function(v) v  # Define the preconditioner function Sig as the identity function
A_func <- function(v) A %*% v  # Define the function A(v) to apply the matrix A to a vector v

bayescg_results <- bayescg(A_func, b, x0, Sig, max_it = 10, tol = 1e-6, reorth = TRUE)

# Find solution using implemented R solver ----
R_solver_results <- solve(A, b)

# Compare results ----
cat('Final solutions for x:\n')

cat('... using Conjugate Gradient Algorithm:\n')
cg_results$x

cat('... using Preconditioned Conjugate Gradient Algorithm (Jacobi):\n')
pcg_results1$x

cat('... using Preconditioned Conjugate Gradient Algorithm (Identity):\n')
pcg_results2$x

cat('... using Bayesian Conjugate Gradient Algorithm:\n')
bayescg_results$x

cat('... using R Solver:\n')
R_solver_results

all.equal(cg_results$x, R_solver_results, tolerance = 1e-5)
all.equal(pcg_results1$x, R_solver_results, tolerance = 1e-5)
all.equal(pcg_results2$x, R_solver_results, tolerance = 1e-5)
all.equal(bayescg_results$x, R_solver_results, tolerance = 1e-5)



