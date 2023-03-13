library(mvtnorm, quietly=T)

# Rcpp packages
library(Rcpp)
library(RcppArmadillo)
library(RcppDist, quietly = T)
sourceCpp("trial1_routine_c.cpp")

# Needed for OpenMP C++ parallel
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")

update_f_j = function(pars, par_index, K, j) {
    
}

update_sigma2 = function(pars, par_index) {
    
}


# -----------------------------------------------------------------------------
# The mcmc routine for samping the parameters
# -----------------------------------------------------------------------------
mcmc_routine = function( y, x, id, init_par, prior_par, par_index,
                         steps, burnin, n_cores, ind){
    
    cl <- makeCluster(n_cores, outfile="")
    registerDoParallel(cl)
    
    pars = init_par
    n = length(y)
    n_par = length(pars)
    chain = matrix( 0, steps, n_par)
    
    group = list(c(par_index$beta), c(par_index$s))
    n_group = length(group)
    
    # proposal covariance and scale parameter for Metropolis step
    pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(group[[j]]))*0.001
    pscale = rep( 1, n_group)
    # load('Model_out/mcmc_out_2_21.rda')
    # pcov = mcmc_out$pcov
    # pscale = mcmc_out$pscale
    # rm(mcmc_out)

    accept = rep( 0, n_group)
    EIDs = unique(id)
    
    # Begin the MCMC algorithm --------------------------------------------------
    chain[1,] = pars
    for(ttt in 2:steps){
        
        # Setting up K_1, K_2
        K = update_K_j(pars, par_index, length(unique(x)))
        K_1 = K[[1]]
        K_2 = K[[2]]
        
        # Gibbs update for f1 and f2
        pars[par_index$f1] = update_f_j(pars, par_index, K, 1)
        pars[par_index$f2] = update_f_j(pars, par_index, K, 2)
        
        
        # Evaluate the log_post of the initial pars
        log_post_prev = fn_log_post_continuous( pars, prior_par, par_index, y, x, id, K)
        
        if(!is.finite(log_post_prev)){
            print("Infinite log-posterior; choose better initial parameters")
            break
        }
        
        for(j in 1:n_group){
            
            # Propose an update
            ind_j = group[[j]]
            proposal = pars
            proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],sigma=pcov[[j]]*pscale[j])
            
            # Compute the log density for the proposal
            log_post = fn_log_post_continuous(proposal, prior_par, 
                                              par_index, y_1, y_2, 
                                              t, id)
            
            # Only propose valid parameters during the burnin period
            if(ttt < burnin){
                while(!is.finite(log_post)){
                    print('bad proposal')
                    proposal = pars
                    proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],
                                               sigma=pcov[[j]]*pscale[j])
                    log_post = fn_log_post_continuous(proposal, prior_par, 
                                                      par_index, y_1, y_2, 
                                                      t, id)
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
        
        if(ttt%%1==0)  cat('--->',ttt,'\n')
    }
    # ---------------------------------------------------------------------------
    
    stopCluster(cl)
    print(accept/(steps-burnin))
    
    return(list( chain=chain[burnin:steps,], accept=accept/(steps-burnin),
                 pscale=pscale, pcov = pcov))
}
# -----------------------------------------------------------------------------
