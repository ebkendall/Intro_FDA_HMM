library(latex2exp)

# -----------------------------------------------------------------------------
# (1) Determining the optimal K for the B-spline approach using an analogue 
#       to adjusted R^2
# -----------------------------------------------------------------------------
load('Data/big_B.rda')
load('Data/data_format_1.rda')
load('Data/y_mat_1.rda')
K = length(big_B)
N = nrow(y_mat)
n_i = length(unique(data_format[,"t"]))

metric_val = rep(0, K)

for(k in 1:K) {
    print(k)
    
    # Taking the results from multiple seeds applied to the same dataset
    chain_list = vector(mode = "list", length = 3)
    for(j in 1:3) {
        load(paste0('Model_out/mcmc_out_', j, '_', k, '_bspline.rda'))
        chain_list[[j]] = mcmc_out$chain[1:5000, ]
        par_index = mcmc_out$par_index
    }
    
    stacked_chains = do.call( rbind, chain_list)
    
    # posterior mean
    beta_1 = colMeans(stacked_chains[,par_index$beta_1])
    beta_2 = colMeans(stacked_chains[,par_index$beta_2])
    
    # posterior mode
    z_1 = as.numeric(colMeans(stacked_chains[,par_index$Z_1]) > 0.5)
    z_2 = as.numeric(colMeans(stacked_chains[,par_index$Z_2]) > 0.5)
    
    xi_1 = z_1 * beta_1
    xi_2 = z_2 * beta_2
    
    sum_val = 0
    for(i in 1:N) {
        # Assuming labels are known
        state_lab = data_format[data_format[,"id"] == i, "true_state"]
        B_1 = B_2 = big_B[[k]]
        B_1[state_lab == 2, ] = 0
        B_2[state_lab == 1, ] = 0
        y_hat = t(B_1 %*% xi_1 + B_2 %*% xi_2)
        numerator = (n_i - 1) * (y_mat[i, ] - y_hat) %*% t(y_mat[i, ] - y_hat)
        
        denominator = (n_i - mean(c(sum(xi_1 != 0), sum(xi_2 != 0)))) * 
                        ((y_mat[i, ,drop=F] - mean(y_mat[i, ])) %*% t(y_mat[i, ,drop=F] - mean(y_mat[i, ])))
        sum_val = sum_val + (numerator / denominator)
    }
    
    metric_val[k] = 1 - (1/N) * sum_val
}

K_values = sapply(big_B, ncol)
print("Metric for determining the best K: ")
print(paste0(K_values, ": ", round(metric_val, digits = 4)))

# -----------------------------------------------------------------------------
# (2) Determining the optimal percentage of variability explained for FPCA 
# -----------------------------------------------------------------------------
P_1 = c(0.99, 0.95, 0.9, 0.85, 0.8, 0.75)
metric_val_fpca = rep(0, length(P_1))

for(k in 1:length(P_1)) {
    print(k)
    
    # Taking the results from multiple seeds applied to the same dataset
    chain_list = vector(mode = "list", length = 3)
    for(j in 1:3) {
        load(paste0('Model_out/mcmc_out_', j, '_', k, '_fpca.rda'))
        chain_list[[j]] = mcmc_out$chain[1:5000, ]
        par_index = mcmc_out$par_index
    }
    stacked_chains = do.call( rbind, chain_list)
    
    # posterior mean
    beta_1 = colMeans(stacked_chains[,par_index$beta_1])
    beta_2 = colMeans(stacked_chains[,par_index$beta_2])
    
    sum_val = 0
    for(i in 1:N) {
        # Assuming labels are known
        state_lab = data_format[data_format[,"id"] == i, "true_state"]
        B_1 = mcmc_out$B_1
        B_2 = mcmc_out$B_2
        B_1[state_lab == 2, ] = 0
        B_2[state_lab == 1, ] = 0
        y_hat = t(B_1 %*% beta_1 + B_2 %*% beta_2)
        numerator = (n_i - 1) * (y_mat[i, ] - y_hat) %*% t(y_mat[i, ] - y_hat)
        
        denominator = (n_i - mean(c(sum(beta_1 != 0), sum(beta_2 != 0)))) * 
            ((y_mat[i, ,drop=F] - mean(y_mat[i, ])) %*% t(y_mat[i, ,drop=F] - mean(y_mat[i, ])))
        sum_val = sum_val + (numerator / denominator)
    }
    
    metric_val_fpca[k] = 1 - (1/N) * sum_val
}

print("Metric for determining the best P1: ")
print(paste0(P_1*100, "%: ", round(metric_val_fpca, digits = 4)))

# -----------------------------------------------------------------------------
# (3) Empirical mean square error (EMSE)
# -----------------------------------------------------------------------------
load('Data/true_fnc_vals.rda')
n_sim = 100
EMSE_bspline = EMSE_fpca = vector(mode = "list", length = 2)
EMSE_bspline[[1]] = EMSE_bspline[[2]] = 
        EMSE_fpca[[1]] = EMSE_fpca[[2]] = matrix(nrow = n_sim, ncol = n_i)

for(i in 1:n_sim) {
    print(i)
    
    # Adaptive B-spline approach (using K = 10 based on test (1)) -------------
    load(paste0('Model_out/mcmc_out_', i, '_4_bspline.rda'))
    par_index = mcmc_out$par_index
    
    # posterior mean
    beta_1 = colMeans(mcmc_out$chain[,par_index$beta_1])
    beta_2 = colMeans(mcmc_out$chain[,par_index$beta_2])
    
    # posterior mode
    z_1 = as.numeric(colMeans(mcmc_out$chain[,par_index$Z_1]) > 0.5)
    z_2 = as.numeric(colMeans(mcmc_out$chain[,par_index$Z_2]) > 0.5)
    
    xi_1 = z_1 * beta_1
    xi_2 = z_2 * beta_2
    
    f_1_hat = c(mcmc_out$B_1 %*% xi_1)
    f_2_hat = c(mcmc_out$B_2 %*% xi_2)
    
    EMSE_bspline[[1]][i,] = (f_1_hat - c(fnc_vals[[1]]))^2
    EMSE_bspline[[2]][i,] = (f_2_hat - c(fnc_vals[[2]]))^2
    
    # FPCA approach (using __% based on the test (2)) -------------------------
    load(paste0('Model_out/mcmc_out_', i, '_1_fpca.rda'))
    par_index = mcmc_out$par_index
    
    # posterior mean
    beta_1 = colMeans(mcmc_out$chain[,par_index$beta_1])
    beta_2 = colMeans(mcmc_out$chain[,par_index$beta_2])
    
    f_1_hat = c(mcmc_out$B_1 %*% beta_1)
    f_2_hat = c(mcmc_out$B_2 %*% beta_2)
    
    EMSE_fpca[[1]][i,] = (f_1_hat - c(fnc_vals[[1]]))^2
    EMSE_fpca[[2]][i,] = (f_2_hat - c(fnc_vals[[2]]))^2
    
}

EMSE_1_bspline = colMeans(EMSE_bspline[[1]])
EMSE_2_bspline = colMeans(EMSE_bspline[[2]])
plot(EMSE_1_bspline, type = 'l', xlab = "t", ylab = "EMSE(t)",
     ylim = c(0, max(EMSE_1_bspline, EMSE_2_bspline)))
lines(EMSE_2_bspline, lty = 2)
legend("top", legend = c(TeX(r'($g_1$)'), TeX(r'($g_2$)')), lty = c(1, 2), ncol = 2)

EMSE_1_fpca = colMeans(EMSE_fpca[[1]])
EMSE_2_fpca = colMeans(EMSE_fpca[[2]])
plot(EMSE_1_fpca, type = 'l', xlab = "t", ylab = "EMSE(t)",
     ylim = c(0, max(EMSE_1_fpca, EMSE_2_fpca)))
lines(EMSE_2_fpca, lty = 2)
legend("top", legend = c(TeX(r'($g_1$)'), TeX(r'($g_2$)')), lty = c(1, 2), ncol = 2)

# -----------------------------------------------------------------------------
# (4) Mean integrated square error
# -----------------------------------------------------------------------------
MISE_1_bspline = 0.01 * sum(EMSE_1_bspline)
MISE_2_bspline = 0.01 * sum(EMSE_2_bspline)

# -----------------------------------------------------------------------------
# (5) Empirical coverage probabilities for function values, t.p.m, and variance 
#     terms and the posterior means boxplots
# -----------------------------------------------------------------------------
# Coverage probabilities ------------------------------------------------------
load('Data/true_pars.rda')
load('Data/true_par_index.rda')
pars_interest = c(pars[c(par_index$init, par_index$omega, par_index$sigma2)],
                  fnc_vals[[1]], fnc_vals[[2]])
cred_set_bspline = cred_set_fpca = rep(list(matrix(ncol=2,nrow=n_sim)), 
                                            length(pars_interest))
post_means_bspline = post_means_fpca = matrix(nrow = n_sim, ncol = length(pars_interest))

for (i in 1:n_sim) {
    print(i)
    
    # Adaptive B-spline approach (using K = 10 based on test (1)) -------------
    load(paste0('Model_out/mcmc_out_', i, '_4_bspline.rda'))
    par_index = mcmc_out$par_index
    B_1 = mcmc_out$B_1
    B_2 = mcmc_out$B_2
    
    est_line1 = est_line2 = matrix(nrow = 5000, ncol = n_i)
    for(w in 1:5000) {
        beta_1 = mcmc_out$chain[w, par_index$beta_1]
        z_1    = mcmc_out$chain[w, par_index$Z_1]
        beta_1_new = beta_1 * z_1
        
        beta_2 = mcmc_out$chain[w, par_index$beta_2]
        z_2    = mcmc_out$chain[w, par_index$Z_2]
        beta_2_new = beta_2 * z_2
        
        est_line1[w,] = c(B_1 %*% beta_1_new)
        est_line2[w,] = c(B_2 %*% beta_2_new)
    }
    
    new_chain = cbind(mcmc_out$chain[1:5000, c(par_index$init, par_index$omega, par_index$sigma2)], 
                      est_line1, est_line2)
    for(j in 1:length(pars_interest)) {
        cred_set_bspline[[j]][i, 1] = round(quantile( new_chain[,j], prob=.025), 4)
        cred_set_bspline[[j]][i, 2] = round(quantile( new_chain[,j], prob=.975), 4)   
    }
    
    post_means_bspline[i, ] = colMeans(new_chain)
    
    # FPCA approach (using __% based on the test (2)) -------------------------
}

cov_df_bspline = rep(NA, length(pars_interest))
for(i in 1:length(pars_interest)) {
    val = pars_interest[i]
    cov_df_bspline[i] = mean(cred_set_bspline[[i]][,1] <= val & 
                                 val <= cred_set_bspline[[i]][,2], na.rm=T)
}

# Boxplots --------------------------------------------------------------------
labels = c('logit initial', 'omega_1', 'omega_2', 'sigma2',
           paste0('f1(', 1:n_i, ')'), paste0('f2(', 1:n_i, ')'))
pdf('Plots/boxplots_bspline.pdf')
par(mfrow=c(3, 3), mar=c(5, 2.5, 2.5, 2.5))  # mar = c(bottom, left, top, right)
for(r in 1:length(labels)){	
    
    truth <- pars_interest[r]
    covg <- cov_df_bspline[r]
    
    plotRange <- c( min( min(post_means_bspline[,r]), truth), 
                    max( max(post_means_bspline[,r]), truth))
    
    boxplot( post_means_bspline[,r], names=c('Posterior means'), ylim=plotRange, ylab=NA, main=labels[r], 
             outline=TRUE, cex.axis=1.25, cex.main=1.5)
    
    abline(h=truth,col='green',lwd=4,lty=5)
    mtext(paste0('.95 coverage = ',toString(round(covg,3))), side=1, line=2.5, cex=1.25)

}
dev.off()

