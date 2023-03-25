library(fda)

set.seed(102)

# -----------------------------------------------------------------------------
# --------------------- Data Generation ---------------------------------------
# -----------------------------------------------------------------------------

# The same evaluation points will be used for all patients
t = seq(0.01,1,by = 0.01) # NEXT: try t_i = runif(100, 0, 1)
sigma2 = 0.1^2
N = 50

# Generate a B-spline of degree 3 (i.e. order 4) with K=10 basis functions for each state
basis_s1 = create.bspline.basis(range(t), nbasis=10, norder=4)

basis_s1_eval_pts = getbasismatrix(t, basis_s1, nderiv=0) # same as eval.basis(t, basis_s1)

# Whichever beta we set to 0 is an indication that it is not necessary
coef_s1 = c(-2, 0,  1.5, 1.5, 0, -1, -0.5, -1, 0, 0) # beta_1

true_fnc_1 = basis_s1_eval_pts %*% coef_s1

# Generate the data -----------------------------------------------------------
y = matrix(nrow = nrow(true_fnc_1), ncol = N)
for(i in 1:N) {
    y[,i] = true_fnc_1 + rnorm(length(c(true_fnc_1)), mean = 0, sd = sqrt(sigma2))
}

# -----------------------------------------------------------------------------
# ---------------------- Data Smoothing ---------------------------------------
# -----------------------------------------------------------------------------

# Create an initial basis based on some number of basis that you want
basis_num = 15
spline_est = create.bspline.basis(range(t), nbasis=basis_num, norder=4)

# Evaluate those basis functions at the time points observed to form the matrix
# (n_i x basis_num)
basismat = eval.basis(t, spline_est)

# Using simple least squares, fit the basis to data to smooth things
# We get the coefficient matrix for each basis for each subject (basis_num x N)
spline_coeff = lsfit(basismat, y, intercept=FALSE)$coef
# spline_coeff = basismat \ y

# Turn the data into estimated functions to then apply FPCA
y_func = smooth.basis(t, y, spline_est)

# Plot to see results
plot(y[,1])
fitted_curve = basismat %*% y_func$fd$coefs[,1]
lines(fitted_curve, col = 'red')

# Saving the functional object
y_func.fd = y_func$fd

# Looking at the estimated covariance
y_func.bifd = var.fd(y_func.fd)
y_func_var = eval.bifd(t, t, y_func.bifd) 

# We can compare these results to doing this formula by hand
y_est = basismat %*% y_func$fd$coefs
y_bars = rowMeans(y_est)

K = matrix(nrow = nrow(y_est), ncol = nrow(y_est))
for(s in 1:nrow(K)) {
    for(t in 1:ncol(K)) {
        temp_sum = 0
        for(i in 1:N) {
            temp_sum = temp_sum + (y_est[s, i] - y_bars[s]) * (y[t, i] - y_bars[t])
        }
        K[s,t] = (1/(N)) * temp_sum
    }
}

# We still need to account for the variance not explained by the smoothing of the data
# This is the residual variance-covariance matrix
# different attempt
y_est = t(y_func$fd$coefs) %*% t(basismat)
K = (1/N) * (t(y_est) %*% y_est)

# -----------------------------------------------------------------------------
# ------------------------------- FPCA ----------------------------------------
# -----------------------------------------------------------------------------

# Simplest Approach from Rao 1958
y = t(y)
V = (1/N) * (t(y) %*% y)
w = 0.01

eigen_V = eigen(V)
rho = eigen_V$values * w
xi = eigen_V$vectors * (1/sqrt(w))

round(rho/sum(rho), digits = 3)

# Tidemann miller
# define the basis needed and then we will do MCMC on the coefficients
psi = cbind(rho[1] * xi[,1], rho[2] * xi[,2])
