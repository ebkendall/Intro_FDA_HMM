# This script file produces trace plots and histograms of the mcmc output files
library(tidyverse)
library(gridExtra)
library(latex2exp)

dir = 'Model_out/'

args = commandArgs(TRUE)

spline_or_fpca = as.numeric(args[1]) 

# Size of posterior sample from mcmc chains
n_post = 5000
# Step number at which the adaptive tuning scheme was frozen
burnin = 5000
# Total number of steps the mcmc algorithm is computed for
steps = 10000
# Matrix row indices for the posterior sample to use
index_post = (steps - burnin - n_post + 1):(steps - burnin)

index_seeds = c(1:3)
trial_num = 4

if(spline_or_fpca == 1) {
    load(paste0('Model_out/mcmc_out_', index_seeds[1], '_', trial_num, '_bspline.rda'))
    l = ncol(mcmc_out$B_1)
    labels = c('logit initial', 'omega_1', 'omega_2', 'sigma2',
            paste0('beta1(', 1:l, ')'), paste0('beta2(', 1:l, ')'),
            paste0('theta1(', 1:l, ')'), paste0('theta2(', 1:l, ')'),
            paste0('Z1(', 1:l, ')'), paste0('Z2(', 1:l, ')'))
} else if(spline_or_fpca == 0) {
    load(paste0('Model_out/mcmc_out_', index_seeds[1], '_', trial_num, '_fpca.rda'))
    l = ncol(mcmc_out$B_1)
    labels = c('logit initial', 'omega_1', 'omega_2', 'sigma2',
            paste0('beta1(', 1:l, ')'), paste0('beta2(', 1:l, ')'))
}

load('Data/true_pars.rda')

# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------
chain_list = vector(mode = "list", length = length(index_seeds))
post_means = matrix(nrow = length(index_seeds), ncol = length(labels))

ind = 0

for(seed in index_seeds){

    if(spline_or_fpca == 1) {
        file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_bspline.rda')
    } else if(spline_or_fpca == 0) {
        file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_fpca.rda')
    }
    
    if (file.exists(file_name)) {
        load(file_name)
        ind = ind + 1

        print(mcmc_out$accept)

        # Thinning the chain
        main_chain = mcmc_out$chain[index_post,]
        ind_keep = seq(1, nrow(main_chain), by=10)

      	chain_list[[ind]] = main_chain[ind_keep, ]
    	post_means[ind,] <- colMeans(main_chain[ind_keep, ])
    }
}

# Plot and save the mcmc trace plots and histograms.
pdf_name=NULL
if(spline_or_fpca == 1) {
    pdf_name = paste0('Plots/mcmc_out_', trial_num, '_bspline.pdf')
} else if(spline_or_fpca == 0) {
    pdf_name = paste0('Plots/mcmc_out_', trial_num, '_fpca.pdf')
}

pdf(pdf_name)
par(mfrow=c(4, 2))

stacked_chains = do.call( rbind, chain_list)
par_mean = par_median = upper = lower = rep( NA, length(labels))

for(r in 1:length(labels)){
    
    plot( NULL, xlab=paste0("true val: ", round(pars[r], 3)), ylab=NA, main=labels[r], xlim=c(1,nrow(chain_list[[1]])),
          ylim=range(stacked_chains[,r]) )
    
    for(seed in 1:length(chain_list)) lines( chain_list[[seed]][,r], type='l', col=seed)

    par_mean[r] = round( mean(stacked_chains[,r]), 4)
    par_median[r] = round( median(stacked_chains[,r]), 4)
    upper[r] = round( quantile( stacked_chains[,r], prob=.975), 4)
    lower[r] = round( quantile( stacked_chains[,r], prob=.025), 4)
    print(paste0(r, ". ", labels[r],": [", lower[r], ", ", upper[r], "]"))

    hist( stacked_chains[,r], breaks=sqrt(nrow(stacked_chains)), ylab=NA, main=NA,
            freq=F, xlab=paste0('Mean = ',toString(par_mean[r]),
                                ' Median = ',toString(par_median[r])))
    abline( v=upper[r], col='red', lwd=2, lty=2)
    abline( v=lower[r], col='purple', lwd=2, lty=2)
    
    abline( v=pars[r], col='green', lwd=2, lty=2)        


}

# Should plot_CI to compare

# means1 = colMeans(stacked_chains[,par_index$beta_1])
# q_95_1 = apply(stacked_chains, 2, quantile, prob = 0.975)[par_index$beta_1]
# q_5_1 = apply(stacked_chains, 2, quantile, prob = 0.05)[par_index$beta_1]
# 
# means2 = colMeans(stacked_chains[,par_index$beta_2])
# q_95_2 = apply(stacked_chains, 2, quantile, prob = 0.975)[par_index$beta_2]
# q_5_2 = apply(stacked_chains, 2, quantile, prob = 0.05)[par_index$beta_2]
# 
# fitted_m1 = B_1 %*% means1
# fitted_q95_1 = B_1 %*% q_95_1
# fitted_q5_1 = B_1 %*% q_5_1
# 
# fitted_m2 = B_2 %*% means2
# fitted_q95_2 = B_2 %*% q_95_2
# fitted_q5_2 = B_2 %*% q_5_2
# 
# points(fitted_m1, col = 'green', pch = 2)
# points(fitted_q95_1, col = 'green', pch = 1)
# points(fitted_q5_1, col = 'green', pch = 1)
# 
# points(fitted_m2, col = 'red', pch = 2)
# points(fitted_q95_2, col = 'red', pch = 1)
# points(fitted_q5_2, col = 'red', pch = 1)

dev.off()
