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

xi_1_hat = xi_2_hat = vector(mode = 'list', length = K)
metric_val = rep(0, K)

for(k in 1:K) {
    print(k)
    load(paste0('Model_out/mcmc_out_1_', k, '_bspline.rda'))
    par_index = mcmc_out$par_index
    
    # posterior mean
    beta_1 = colMeans(mcmc_out$chain[,par_index$beta_1])
    beta_2 = colMeans(mcmc_out$chain[,par_index$beta_2])
    
    # posterior mode
    z_1 = as.numeric(colMeans(mcmc_out$chain[,par_index$Z_1]) > 0.5)
    z_2 = as.numeric(colMeans(mcmc_out$chain[,par_index$Z_2]) > 0.5)
    
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

# -----------------------------------------------------------------------------
# (2) Determining the optimal percentage of variability explained for FPCA 
# -----------------------------------------------------------------------------

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
    
}

EMSE_1_bspline = colMeans(EMSE_bspline[[1]])
EMSE_2_bspline = colMeans(EMSE_bspline[[2]])
plot(EMSE_1_bspline, type = 'l', xlab = "t", ylab = "EMSE(t)",
     ylim = c(0, max(EMSE_1_bspline, EMSE_2_bspline)))
lines(EMSE_2_bspline, lty = 2)
legend("top", legend = c(TeX(r'($g_1$)'), TeX(r'($g_2$)')), lty = c(1, 2), ncol = 2)

# -----------------------------------------------------------------------------
# Mean integrated square error
# -----------------------------------------------------------------------------
MISE_1_bspline = 0.01 * sum(EMSE_1_bspline)
MISE_2_bspline = 0.01 * sum(EMSE_2_bspline)

# -----------------------------------------------------------------------------
# Empirical coverage probabilities for function values, t.p.m, and variance terms
#     and the posterior means boxplots
# -----------------------------------------------------------------------------


