library(mvtnorm, quietly=T)

# Rcpp packages
library(Rcpp)
library(RcppArmadillo)
library(RcppDist, quietly = T)
sourceCpp("mcmc_routine_c.cpp")

# Needed for OpenMP C++ parallel
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")


# -----------------------------------------------------------------------------
# The mcmc routine for samping the parameters
# -----------------------------------------------------------------------------
mcmc_routine = function( y, t, id, init_par, prior_par, par_index,
                         steps, burnin, ind, init_state, B_1_master, B_2_master){
    
    pars = init_par
    n = length(y)
    n_par = length(pars)
    chain = matrix( 0, steps, n_par)
    B_chain = matrix( 0, steps - burnin, length(y))
    
    group = list(c(par_index$init, par_index$omega))
    n_group = length(group)
    
    # proposal covariance and scale parameter for Metropolis step
    pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(group[[j]]))*0.001
    pscale = rep( 1, n_group)
    
    accept = rep( 0, n_group)
    EIDs = unique(id)
    
    # Initialize and store proposed state sequences
    B = list()
    for(i in 1:length(EIDs)) {
        state_sub = init_state[id == EIDs[i]]
        b_temp = matrix(state_sub, ncol = 1)
        B[[i]] = b_temp
    }
    
    # Initialize W_1 and W_2
    W_1 = W_2 = vector(mode = 'list', length = length(EIDs))
    for(i in 1:length(EIDs)) {
        W_1[[i]] = as.numeric(B[[i]] == 1)
        W_2[[i]] = as.numeric(B[[i]] == 2)
    }
    
    # Initialize the B_1 and B_2
    both_B = update_both_B(pars, par_index, B_1_master, B_2_master)
    B_1 = both_B[[1]]
    B_2 = both_B[[2]]
    
    # Begin the MCMC algorithm --------------------------------------------------
    chain[1,] = pars
    for(ttt in 2:steps){
        
        # Gibbs update for beta_1 and beta_2
        pars[par_index$beta_1] = update_beta_j(pars, par_index, W_1, W_2, 
                                               y, id, 1, EIDs, B_1, B_2)
        pars[par_index$beta_2] = update_beta_j(pars, par_index, W_1, W_2, 
                                               y, id, 2, EIDs, B_1, B_2)
        
        # Gibbs update for sigma2
        pars[par_index$sigma2] = update_sigma2(pars, par_index, id, y, EIDs,
                                               W_1, W_2, B_1, B_2)
        
        # B_chain: Metropolis-within-Gibbs update
        B = update_b_i_cpp(pars, par_index, y, id, EIDs, B, B_1, B_2)
        
        # Refresh W
        for(i in 1:length(EIDs)) {
            W_1[[i]] = as.numeric(B[[i]] == 1)
            W_2[[i]] = as.numeric(B[[i]] == 2)
        }
        
        # Update Z_1 and Z_2
        pars[par_index$Z_1] = update_Z_j(pars, prior_par, par_index, y, id, EIDs, B_1_master, B_2_master, 1)
        pars[par_index$Z_2] = update_Z_j(pars, prior_par, par_index, y, id, EIDs, B_1_master, B_2_master, 2)
        
        # Refresh what B_1 and B_2 are after updating Z_1 and Z_2
        both_B = update_both_B(pars, par_index, B_1_master, B_2_master)
        B_1 = both_B[[1]]
        B_2 = both_B[[2]]
        
        # Update theta_1 and theta_2
        # pars[par_index$theta_1] = update_theta_j(pars, par_index, 1)
        # pars[par_index$theta_2] = update_theta_j(pars, par_index, 2)
        
        # Evaluate the log_post of the initial pars
        log_post_prev = fn_log_post_continuous( pars, prior_par, par_index, y, id, EIDs, B_1, B_2)
        
        if(!is.finite(log_post_prev)){
            print("Infinite log-posterior")
            break
        }
        
        chain[ttt, ] = pars
        for(j in 1:n_group){
            
            # Propose an update
            ind_j = group[[j]]
            proposal = pars
            proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],sigma=pcov[[j]]*pscale[j])
            
            # Compute the log density for the proposal
            log_post = fn_log_post_continuous(proposal, prior_par, par_index, y, id, EIDs, B_1, B_2)
            
            # Only propose valid parameters during the burnin period
            if(ttt < burnin){
                while(!is.finite(log_post)){
                    print('bad proposal')
                    proposal = pars
                    proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],
                                               sigma=pcov[[j]]*pscale[j])
                    log_post = fn_log_post_continuous(proposal, prior_par, par_index, y, id, EIDs, B_1, B_2)
                }
            }
            
            # Evaluate the Metropolis-Hastings ratio
            if( log_post - log_post_prev > log(runif(1,0,1)) ){
                log_post_prev = log_post
                pars[ind_j] = proposal[ind_j]
                accept[j] = accept[j] +1
            }
            chain[ttt,ind_j] = pars[ind_j]
            
            # Proposal tuning scheme ------------------------------------------------
            if(ttt < burnin){
                # During the burnin period, update the proposal covariance in each step
                # to capture the relationships within the parameters vectors for each
                # transition.  This helps with mixing.
                if(ttt == 100)  pscale[j] = 1
                
                if(100 <= ttt & ttt <= 2000){
                    temp_chain = chain[1:ttt,ind_j]
                    pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
                    
                } else if(2000 < ttt){
                    temp_chain = chain[(ttt-2000):ttt,ind_j]
                    pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
                }
                if( sum( is.na(pcov[[j]]) ) > 0)  pcov[[j]] = diag( length(ind_j) )
                
                # Tune the proposal covariance for each transition to achieve
                # reasonable acceptance ratios.
                if(ttt %% 30 == 0){
                    if(ttt %% 480 == 0){
                        accept[j] = 0
                        
                    } else if( accept[j] / (ttt %% 480) < .4 ){
                        pscale[j] = (.75^2)*pscale[j]
                        
                    } else if( accept[j] / (ttt %% 480) > .5 ){
                        pscale[j] = (1.25^2)*pscale[j]
                    }
                }
            }
            # -----------------------------------------------------------------------
        }
        # Restart the acceptance ratio at burnin.
        if(ttt == burnin)  accept = rep( 0, n_group)
        if(ttt > burnin) {
            B_chain[ttt - burnin, ] = do.call( 'c', B)
        }
        
        if(ttt%%1==0)  cat('--->',ttt,'\n')
    }
    # ---------------------------------------------------------------------------
    
    print(accept/(steps-burnin))
    
    return(list( chain=chain[burnin:steps,], accept=accept/(steps-burnin),
                 pscale=pscale, pcov = pcov, B_chain = B_chain,
                 B_1 = B_1_master, B_2 = B_2_master, par_index = par_index))
}
# -----------------------------------------------------------------------------

