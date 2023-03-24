library(mvtnorm, quietly=T)

# Rcpp packages
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
sourceCpp("final_routine_c.cpp")

# Needed for OpenMP C++ parallel
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")

# R functions ------------------------------------------------------------------
update_K_j_r = function(pars, par_index, x) {
    s1 = pars[par_index$s][1]
    s2 = pars[par_index$s][2]
    
    U_1 = 1/(sqrt(2 * pi * (s1^2)))
    U_2 = 1/(sqrt(2 * pi * (s2^2)))
    
    # covariance for f_1 and f_2 
    # U_j exp(-(x-t)^2  /  2s_j^2)
    K_1 = K_2 = matrix(nrow=length(x), ncol = length(x))
    
    for(i in 1:length(x)) {
        for(j in 1:length(x)) {
            l = x[i]
            m = x[j]
            K_1[i, j] = U_1 * exp(-((l - m)^2) / (2 * s1^2))
            K_2[i, j] = U_2 * exp(-((l - m)^2) / (2 * s2^2))
        }
    }
    
    return(list(K_1 = K_1,
                K_2 = K_2))
}

fn_log_post_continuous_r <- function(pars, prior_par, par_index, y, x, id, K, EIDs) {
    
    init_logit = c( 1, exp(pars[par_index$init][1]))
    
    # Initial state probabilities
    init = init_logit / sum(init_logit)
    
    # Transition probability matrix
    P = matrix(c(1, exp(pars[par_index$beta][1]),
                   exp(pars[par_index$beta][2]), 1), nrow=2, byrow = T)
    P = P / rowSums(P)
    
    # Parallelized computation of the log-likelihood
    log_total_val = 0
    for(i in unique(id)) {
                                
        f_i = val = 1
        y_i = y[id == i, ,drop = F] # value of response
        
        d_1 = dnorm(x = y_i[1,], mean = pars[par_index$f1][1], sd = sqrt(pars[par_index$sigma2]))
        d_2 = dnorm(x = y_i[1,], mean = pars[par_index$f2][1], sd = sqrt(pars[par_index$sigma2]))
        
        f_i = init %*% diag(c(d_1, d_2))
        
        log_norm = 0
        
        for(k in 2:length(y_i)) {
            
            d_1 = dnorm(x = y_i[k,], mean = pars[par_index$f1][k], sd = sqrt(pars[par_index$sigma2]))
            d_2 = dnorm(x = y_i[k,], mean = pars[par_index$f2][k], sd = sqrt(pars[par_index$sigma2]))
            
            D_i = diag(c(d_1, d_2))
            
            f_i = f_i %*% P %*% D_i
        }
        
        log_total_val = log_total_val + log(sum(f_i))
    }
    
    # Priors for Metropolis Step
    mean = prior_par$prior_mean
    sd = diag(prior_par$prior_sd)
    log_prior_dens = dmvnorm( x=c(pars[par_index$init], pars[par_index$beta], pars[par_index$s]),
                              mean=mean, sigma=sd, log=T)
    
    # Likelihood contribution from the rest
    log_f1_dens = dmvnorm(x = pars[par_index$f1], mean = rep(0, length(par_index$f1)), sigma = K[[1]], log = T)
    log_f2_dens = dmvnorm(x = pars[par_index$f2], mean = rep(0, length(par_index$f2)), sigma = K[[2]], log = T)

    return(log_total_val + log_prior_dens + log_f1_dens + log_f2_dens)
    
}
# ------------------------------------------------------------------------------

# loading data -----------------------------------------------------------------
load('Data/data_format.rda')
load('Data/true_pars.rda')
load('Data/par_index.rda')
# ------------------------------------------------------------------------------

prior_mean = rep(0, 5)
prior_sd = rep(5, 5)
prior_par = data.frame( prior_mean= prior_mean,
                        prior_sd= prior_sd)

init_state = data_format[,"true_state"]

temp_data = as.matrix(data_format); rownames(temp_data) = NULL
id = temp_data[,"id", drop=F]
y = temp_data[,"y", drop = F]
EIDs = unique(id)

B = list()
for(i in 1:length(EIDs)) {
    state_sub = init_state[id == EIDs[i]]
    b_temp = matrix(state_sub, ncol = 1)
    B[[i]] = b_temp
}


test_functions(pars, par_index)