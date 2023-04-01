library(matrixStats)
library(plotrix)

args <- commandArgs(TRUE)
set.seed(args[1])

spline_or_fpca = as.numeric(args[2])

trial_num = 1

# Loading data
load('Data/data_format.rda')
load('Data/par_index.rda')
load('Data/true_pars.rda')
load('Data/big_B.rda')
load('Data/fnc_vals.rda')

if(spline_or_fpca == 1) {
    load(paste0('Model_out/mcmc_out_',toString(args[1]), '_', trial_num, '_bspline.rda'))
} else if(spline_or_fpca == 0) {
    load(paste0('Model_out/mcmc_out_',toString(args[1]), '_', trial_num, '_fpca.rda'))
}

mcmc_out$B_chain = mcmc_out$B_chain[1:5000, ]
# f1_chain = mcmc_out$chain[15000:25000, par_index$f1]
# f2_chain = mcmc_out$chain[15000:25000, par_index$f2]
true_line1 = fnc_vals$true_fnc_1
true_line2 = fnc_vals$true_fnc_2

# Get estimates of the functions
t_pts = 100
est_line1 = est_line2 = matrix(nrow = 5000, ncol = t_pts)
if(spline_or_fpca == 1) {
    B_1 = big_B[['10']]
    B_2 = big_B[['10']]
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
        beta_1 = mcmc_out$chain[i, 5:23]
        beta_2 = mcmc_out$chain[i, 24:42]

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
    
    
    # f1 --------------------------------------------------------------
    # f1_upper = colQuantiles( f1_chain, probs=.975)
    # f1_lower = colQuantiles( f1_chain, probs=.025)
    # f2_upper = colQuantiles( f2_chain, probs=.975)
    # f2_lower = colQuantiles( f2_chain, probs=.025)
    # 
    # ylimit = c(min(min(y_i), min(f1_lower, f2_lower)), max(max(y_i), max(f1_upper, f2_upper)))
    # 
    # plotCI( x=colMeans(f1_chain), ui=f1_upper, li=f1_lower, ylim = ylimit,
    #         main=paste0('Function 1 (green) and 2 (purple)'), xlab='time', ylab=NA, xaxt='n', col='turquoise',
    #         col.axis='green', pch=20, cex=1, sfrac=.0025)
    # 
    # plotCI( x=colMeans(f2_chain), ui=f2_upper, li=f2_lower, ylim = ylimit, col='pink',
    #         pch=20, cex=1, sfrac=.0025, add = T)
    # 
    # lines(pars[par_index$f1], col = "blue")
    # lines(pars[par_index$f2], col = "red")
    # points(y_i, col = col_i, pch = 17, cex = 2)
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
