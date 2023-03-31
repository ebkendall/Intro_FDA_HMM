source("spline_routine.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

load('Data/big_B.rda')
# load(paste0('Data/data_format_', ind, '.rda'))
load(paste0('Data/data_format_1.rda'))

for(k in 1:length(big_B)) {
    
    # These are initialized such that all Z_1 = Z_2 = 1
    B_1_master = B_2_master = big_B[[k]]
    
    # Number of B-splines
    K = ncol(big_B[[k]])
    
    print(paste0("Number of Splines: ", K))
    
    # Setting the initial parameter values
    init_par = c(     0,                   # init
                 -4, -4,                   # omega
                   0.01,                   # sigma2
                 rep(0  , K), rep(0  , K), # beta_1 & beta_2
                 rep(0.5, K), rep(0.5, K), # theta_1 & theta_2
                 rep(1  , K), rep(1  , K)) # Z_1 & Z_2
    
    # Setting the par_index for easy access in the algorithm
    par_index = list( init=1, omega = 2:3, sigma2 = 4, 
                      beta_1 = 5:(5 + K - 1), beta_2 = (5 + K):(5 + 2*K - 1),
                      theta_1 = (5 + 2*K):(5 + 3*K - 1), theta_2 = (5 + 3*K):(5 + 4*K - 1), 
                      Z_1 = (5 + 4*K):(5 + 5*K - 1), Z_2 = (5 + 5*K):(5 + 6*K - 1))
    
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
                            steps, burnin, ind, init_state, B_1_master, B_2_master)
    e_time = Sys.time() - s_time; print(e_time)
    
    save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", k, "_bspline.rda"))
}