# This script file produces trace plots and histograms of the mcmc output files
library(tidyverse)
library(gridExtra)
library(latex2exp)

dir = 'Model_out/' 
args <- commandArgs(TRUE)

# Size of posterior sample from mcmc chains
n_post = 5000
# Step number at which the adaptive tuning scheme was frozen
burnin = 5000
# Total number of steps the mcmc algorithm is computed for
steps = 30000
# Matrix row indices for the posterior sample to use
index_post = (steps - burnin - n_post + 1):(steps - burnin)

index_seeds = c(1:3)

labels = c('logit initial', 'beta_1', 'beta_2', 'sigma2', 's1', 's2',
           paste0('f1(x', 1:10, ')'), paste0('f2(x', 1:10, ')'))
load('Data/par_index.rda')
load('Data/true_pars.rda')

# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------
chain_list = vector(mode = "list", length = length(index_seeds))
post_means = matrix(nrow = length(index_seeds), ncol = length(labels))

ind = 0

for(seed in index_seeds){

    file_name = paste0(dir,'mcmc_out_',toString(seed), '.rda')
    
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

pdf(paste0('Plots/mcmc_out_', 1, '.pdf'))
par(mfrow=c(4, 2))

stacked_chains = do.call( rbind, chain_list)
par_mean = par_median = upper = lower = rep( NA, length(labels))

labels_sub <- 1:length(labels)

for(r in 1:length(labels_sub)){

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

dev.off()
