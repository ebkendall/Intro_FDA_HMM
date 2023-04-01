source("fpca_routine.r")

library(expm)

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

converge_or_sim = as.numeric(args[2]) 

if(converge_or_sim == 1) {
    # We run the MCMC for the same data for multiple seeds to check convergence
    # and determine the optimal P_1
    load('Data/data_format_1.rda')
    load('Data/y_mat_1.rda')
    
    # Different percentages to fine tune the PCA methods
    P_1 = c(0.99, 0.95, 0.9, 0.85, 0.8, 0.75)
} else if (converge_or_sim == 0) {
    # We have established the optimal P_1, and now will run for 100 datasets
    load(paste0('Data/data_format_', ind, '.rda'))
    load(paste0('Data/y_mat_', ind, '.rda'))
    P_1 = 0.99
}

# Eigendecomposition for FPCA ----------------------------------------
N = nrow(y_mat)
w = 0.01 # the separation between each time point

y_means = colMeans(y_mat)
y_means = matrix(rep(y_means, N), nrow = N, byrow = T)

V_mat = (1/N) * (t(y_mat - y_means) %*% (y_mat - y_means))
eigen_V = eigen(V_mat)

for(k in 1:length(P_1)) {
    
    print(paste0("Percentage of variation explained: ", P_1[k]))
    
    # Determine the number of eigenfunctions
    rho_k = rep(0, length(eigen_V$values))
    lambda_k = rep(0, length(eigen_V$values))
    for(i in 1:length(rho_k)){
        rho_k[i] = sum(eigen_V$values[1:i]) / sum(eigen_V$values)
        if(eigen_V$values[i] < w) lambda_k[i] = 1
    }
    
    # Determining the number of principal components
    rho_k = as.numeric(rho_k >= P_1[k])
    K = min(min(which(rho_k == 1)), min(which(lambda_k == 1)))
    
    p_comp_keep = 1:K
    
    B_1 = B_2 = eigen_V$vectors[,p_comp_keep]
    
    # Setting the initial parameter values
    init_par = c( 0,         # init
                  -2, -2,    # omega
                  0.01,      # sigma2
                  rep(0, K), # beta_1
                  rep(0, K)) # beta_2
    
    # Initializing the proposed state sequence to the true state sequence
    par_index = list( init=1, omega = 2:3, sigma2 = 4, beta_1 = 5:(5 + K - 1), 
                      beta_2 = (5 + K):length(init_par))
    
    # Initializing the proposed state sequence to the true state sequence
    init_state = data_format[,"true_state"]
    
    # Setting uninformative priors for the Metropolis Step to update init & omega
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
    
    s_time = Sys.time()
    mcmc_out = mcmc_routine(y, t, id, init_par, prior_par, par_index,
                            steps, burnin, ind, init_state, B_1, B_2)
    e_time = Sys.time() - s_time; print(e_time)

    save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", k, "_fpca.rda"))
}



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

