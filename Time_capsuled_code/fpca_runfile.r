source("fpca_routine.r")

library(expm)

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

trial_num = 1

load('Data/data_format.rda')
load('Data/par_index.rda')
load('Data/y_mat.rda')

init_state = data_format[,"true_state"]

prior_mean = rep(0, 3)
prior_sd = rep(20, 3)
prior_par = data.frame( prior_mean= prior_mean,
                        prior_sd= prior_sd)

temp_data = as.matrix(data_format); rownames(temp_data) = NULL
id = temp_data[,"id"]
y = temp_data[,"y"]
t = temp_data[,"t"]
steps = 10000
burnin = 5000

# Establishing the basis from FPCA --------------------------------------------
# B_1 = B_2 = big_B[['10']]
N = nrow(y_mat)
w = diff(t)[1]

y_means = colMeans(y_mat)
y_means = matrix(rep(y_means, 50), nrow = 50, byrow = T)

V_mat = (1/N) * (t(y_mat - y_means) %*% (y_mat - y_means))
eigen_V = eigen(V_mat)

# Determine the number of eigenfunctions
rho_k = rep(0, length(eigen_V$values))
lambda_k = rep(0, length(eigen_V$values))
for(i in 1:length(rho_k)){
    rho_k[i] = sum(eigen_V$values[1:i]) / sum(eigen_V$values)
    if(eigen_V$values[i] < w) lambda_k[i] = 1
}
rho_k = as.numeric(rho_k >= 0.9)
K = min(min(which(rho_k == 1)), min(which(lambda_k == 1)))

p_comp_keep = 1:K

B_1 = B_2 = eigen_V$vectors[,p_comp_keep]

# rho = eigen_V$values * w
# xi = eigen_V$vectors * (1/sqrt(w))
# -----------------------------------------------------------------------------

# # (NEW) Establishing the basis from FPCA --------------------------------------
# # First smooth the observations by fitting a smooth spline function to each observation
# basis_num = 15
# spline_est = create.bspline.basis(range(t), nbasis=basis_num, norder=4)
# 
# # Evaluate those basis functions at the time points observed to form the matrix
# # (n_i x basis_num)
# basismat = eval.basis(unique(t), spline_est)
# basismat = t(basismat)
# 
# # Using least squares to estimate the coefficients "Regression Spline" (FDA w/ R & Matlab)
# spline_coeff = lsfit(basismat, t(y_mat), intercept=FALSE)$coef
# C = t(spline_coeff)
# 
# Y_mat = C %*% basismat
# 
# # Next, W is the numerical integral of the cross product of the B-splines
# # We can approximate this with bins
# W_est = w * (basismat %*% t(basismat))
# W_sqrt = sqrtm(W_est)
# 
# # Then we want the spline estimates for the eigenfunctions
# eigen_mat = (1/N) * (W_sqrt %*% t(C) %*% C %*% W_sqrt)
# eigen_results = eigen(eigen_mat)
# lambda = eigen_results$values
# u = eigen_results$vectors
# b = solve(W_sqrt) %*% u
# -----------------------------------------------------------------------------

init_par = c( 0,      # init
              -2, -2,  # omega
              0.1^2,  # sigma2
              rep(0, length(p_comp_keep)), # beta_1
              rep(0, length(p_comp_keep))) # beta_2

par_index = list( init=1, omega = 2:3, sigma2 = 4, beta_1 = 5:(5+length(p_comp_keep)-1), 
                  beta_2 = (5+length(p_comp_keep)):length(init_par))

s_time = Sys.time()

mcmc_out = mcmc_routine(y, t, id, init_par, prior_par, par_index,
                        steps, burnin, ind, init_state, B_1, B_2)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, "_fpca.rda"))


