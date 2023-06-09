library(matrixStats)
library(plotrix)
args <- commandArgs(TRUE)
set.seed(args[1])

# Loading data
load('Data/data_format.rda')
load('Data/par_index.rda')
load('Data/true_pars.rda')
load('Data/big_B.rda')
load(paste0('Model_out/mcmc_out_',toString(args[1]),'.rda'))

mcmc_out$B_chain = mcmc_out$B_chain[1:5000, ]
# f1_chain = mcmc_out$chain[15000:25000, par_index$f1]
# f2_chain = mcmc_out$chain[15000:25000, par_index$f2]

B1 = B2 = big_B[['10']]
line1 = c(B1 %*% pars[par_index$beta_1])
line2 = c(B2 %*% pars[par_index$beta_2])

simulation=T

EIDs = unique(data_format[,'id'])

# New patients ---------------------------------------------------------------
pdf(paste0('Plots/chart_plot_seed',toString(args[1]), '.pdf'))
panels = c(4, 1)
par(mfrow=panels, mar=c(2,4,2,4), bg='black', fg='green')
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
    plot(y_i, col = col_i, pch = 17, cex = 2, main = "Function 1 (green) and 2 (purple)", xlab = 'time', ylab = 'y')
    lines(line1, col = 4)
    lines(line2, col = 2)
    
    # probs --------------------------------------------------------------
    barplot( rbind( colMeans(mcmc_out$B_chain[, indices_i] == 1),
                    colMeans(mcmc_out$B_chain[, indices_i] == 2)), 
             col=c( 'dodgerblue', 'firebrick1'), 
             xlab='time', xaxt='n', space=0, 
             col.main='green', border=NA) 
    grid( nx=NA, NULL, col='white')
    legend( 'topright', inset=c(0,-.28), xpd=T, horiz=T, bty='n', x.intersp=.75,
            legend=c( 'f1', 'f2'), pch=15, pt.cex=1.5, 
            col=c( 'dodgerblue', 'firebrick1'))
    axis( side=2, at=0:1, col.axis='green')
    
}
dev.off()
