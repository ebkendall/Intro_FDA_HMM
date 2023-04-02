library(matrixStats)
library(plotrix)
library(latex2exp)

args <- commandArgs(TRUE)
set.seed(args[1])

spline_or_fpca = as.numeric(args[2])

# Loading data
load(paste0('Data/data_format_', toString(args[1]), '.rda'))
load('Data/true_pars.rda')
load('Data/big_B.rda')
load('Data/true_fnc_vals.rda')

if(spline_or_fpca == 1) {
    load(paste0('Model_out/Simulation/mcmc_out_',toString(args[1]), '_4_bspline.rda'))
    par_index = mcmc_out$par_index
} else if(spline_or_fpca == 0) {
    load(paste0('Model_out/Simulation/mcmc_out_',toString(args[1]), '_1_fpca.rda'))
    par_index = mcmc_out$par_index
}

mcmc_out$B_chain = mcmc_out$B_chain[1:5000, ]
true_line1 = fnc_vals$true_fnc_1
true_line2 = fnc_vals$true_fnc_2

# Get estimates of the functions
t_pts = 100
est_line1 = est_line2 = matrix(nrow = 5000, ncol = t_pts)
if(spline_or_fpca == 1) {
    B_1 = mcmc_out$B_1
    B_2 = mcmc_out$B_2
    for(i in 1:5000) {
        beta_1 = mcmc_out$chain[i, par_index$beta_1]
        z_1    = mcmc_out$chain[i, par_index$Z_1]
        beta_1_new = beta_1 * z_1

        beta_2 = mcmc_out$chain[i, par_index$beta_2]
        z_2    = mcmc_out$chain[i, par_index$Z_2]
        beta_2_new = beta_2 * z_2

        est_line1[i,] = c(B_1 %*% beta_1_new)
        est_line2[i,] = c(B_2 %*% beta_2_new)
    }
} else if(spline_or_fpca == 0) {
    B_1 = mcmc_out$B_1
    B_2 = mcmc_out$B_2
    for(i in 1:5000) {
        beta_1 = mcmc_out$chain[i, par_index$beta_1]
        beta_2 = mcmc_out$chain[i, par_index$beta_2]

        est_line1[i,] = c(B_1 %*% beta_1)
        est_line2[i,] = c(B_2 %*% beta_2)
    }
}

simulation=T

EIDs = unique(data_format[,'id'])

# New patients ---------------------------------------------------------------
pdf_name=NULL
if(spline_or_fpca == 1) {
    pdf_name = paste0('Plots/chart_plot_seed',toString(args[1]), '_bspline.pdf')
} else if(spline_or_fpca == 0) {
    pdf_name = paste0('Plots/chart_plot_seed',toString(args[1]), '_fpca.pdf')
}

pdf(pdf_name)
panels = c(3, 1)
par(mfrow=panels, mar=c(2,4,2,4))
for(i in EIDs){
    print(which(EIDs == i))
    indices_i = (data_format[,'id']==i)
    
    y_i = data_format[indices_i, "y"]
    t_i = data_format[indices_i, "t"]
    col_i = data_format[indices_i, "true_state"]
    col_i[col_i == 1] = 4
    
    plot(y_i, col = col_i, pch = 17, cex = 2, main = "Function 1 (blue) and 2 (red)", xlab = 'time', ylab = 'y')
    lines(true_line1, col = 4)
    lines(true_line2, col = 2)
    
    # probs --------------------------------------------------------------
    barplot( rbind( colMeans(mcmc_out$B_chain[, indices_i] == 1),
                    colMeans(mcmc_out$B_chain[, indices_i] == 2)), 
             col=c( 'dodgerblue', 'firebrick1'), 
             xlab='time', xaxt='n', space=0, border=NA) 
    grid( nx=NA, NULL, col='white')
    legend( 'topright', inset=c(0,-.28), xpd=T, horiz=T, bty='n', x.intersp=.75,
            legend=c( 'f1', 'f2'), pch=15, pt.cex=1.5, 
            col=c( 4, 2))

    # What is the estimated function -----------------------------------------
    l1_upper = colQuantiles( est_line1, probs=.975)
    l1_lower = colQuantiles( est_line1, probs=.025)
    l2_upper = colQuantiles( est_line2, probs=.975)
    l2_lower = colQuantiles( est_line2, probs=.025)

    ylim = c(min(c(est_line1, est_line2)), max(c(est_line1, est_line2)))
    plotCI( x=colMeans(est_line1), ui=l1_upper, li=l1_lower,
            main=paste0('Estimated functions'), xlab='time', ylab=NA, xaxt='n', col='blue',
                     pch=20, cex=1, sfrac=.0025, ylim = ylim)
    plotCI( x=colMeans(est_line2), ui=l2_upper, li=l2_lower, col='red',
             pch=20, cex=1, sfrac=.0025, add = T)
    lines(true_line1, col = 4, lwd = 0.5)
    lines(true_line2, col = 2, lwd = 0.5)
}
dev.off()


# Drawing specific patient 42
indices_i = (data_format[,'id']==42)

y_i = data_format[indices_i, "y"]
t_i = data_format[indices_i, "t"]
col_i = data_format[indices_i, "true_state"]
col_i[col_i == 1] = 4

ylim = c(min(c(fnc_vals[[1]], fnc_vals[[2]], y_i)), max(c(fnc_vals[[1]], fnc_vals[[2]], y_i)))
png("Plots/sub_42.png", width = 1000, height = 300)
plot(seq(0.01,1,by=0.01), y_i, col = col_i, pch = col_i+13, cex = 2, ylim = ylim,
     main = "Simulated time series", xlab = 'time', ylab = 'y', cex.main = 2)
lines(seq(0.01,1,by=0.01), fnc_vals[[1]], col = 4, lwd = 2)
lines(seq(0.01,1,by=0.01), fnc_vals[[2]], col = 2, lwd = 2)
legend( 'topright', inset=c(.8,-.18), xpd=T, horiz=T, bty='n', x.intersp=.75,
        legend=c( 'State 1', 'State 2'), pch=15, cex=1.5, 
        col=c( 4, 2))
dev.off()
# probs --------------------------------------------------------------
load(paste0('Model_out/Simulation/mcmc_out_',toString(args[1]), '_4_bspline.rda'))
par_index = mcmc_out$par_index
mcmc_out$B_chain = mcmc_out$B_chain[1:5000, ]
png("Plots/sub_42_pp_bspline.png", width = 1000, height = 300)
barplot( rbind( colMeans(mcmc_out$B_chain[, indices_i] == 1),
                colMeans(mcmc_out$B_chain[, indices_i] == 2)), 
         col=c( 'dodgerblue', 'firebrick1'), 
         xlab='time', xaxt='n', space=0, border=NA, main = "Posterior probabilities (aBS)",
         cex.main = 2) 
grid( nx=NA, NULL, col='white')
axis(1, at = 1:100, labels = seq(0.01,1,by=0.01), cex = 1.5)
legend( 'topright', inset=c(.8,-.18), xpd=T, horiz=T, bty='n', x.intersp=.75,
        legend=c( 'State 1', 'State 2'), pch=15, cex=1.5, 
        col=c( 4, 2))
dev.off()

load(paste0('Model_out/Simulation/mcmc_out_',toString(args[1]), '_1_fpca.rda'))
par_index = mcmc_out$par_index
mcmc_out$B_chain = mcmc_out$B_chain[1:5000, ]
png("Plots/sub_42_pp_fpca.png", width = 1000, height = 300)
barplot( rbind( colMeans(mcmc_out$B_chain[, indices_i] == 1),
                colMeans(mcmc_out$B_chain[, indices_i] == 2)), 
         col=c( 'dodgerblue', 'firebrick1'), 
         xlab='time', xaxt='n', space=0, border=NA, main = "Posterior probabilities (FPCA)",
         cex.main = 2) 
grid( nx=NA, NULL, col='white')
axis(1, at = 1:100, labels = seq(0.01,1,by=0.01), cex = 1.5)
legend( 'topright', inset=c(.8,-.18), xpd=T, horiz=T, bty='n', x.intersp=.75,
        legend=c( 'State 1', 'State 2'), pch=15, cex=1.5, 
        col=c( 4, 2))
dev.off()
