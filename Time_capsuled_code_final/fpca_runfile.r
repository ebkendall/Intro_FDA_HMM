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
    print("Running for multiple P_1")
    
    load('Data/data_format_1.rda')
    load('Data/y_mat_smooth_1.rda')
    save_dir = "Model_out/One_data_set/"
    
    # Different percentages to fine tune the PCA methods
    P_1 = c(0.99, 0.95, 0.9, 0.85, 0.8, 0.75)
} else if (converge_or_sim == 0) {
    # We have established the optimal P_1 and now will run for 100 datasets
    print(paste0("Running for single P_1: ", 0.99))
    
    load(paste0('Data/data_format_', ind, '.rda'))
    load(paste0('Data/y_mat_smooth_', ind, '.rda'))
    
    save_dir = "Model_out/Simulation/"
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

    save(mcmc_out, file = paste0(save_dir, "mcmc_out_", ind, "_", k, "_fpca.rda"))
}
