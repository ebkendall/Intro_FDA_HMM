# This workflow file will run everything needed to reproduce the results in 
# the simulation study. 
# Much of this workflow.sh is done in an embarrassingly parallelizable fashion. 

####################
# First generate the data
Rscript data_sim.r

# This code will be used to check convergence as well as determining the 
# optimal K for spline and percentage for FPCA.
# The for-loop above can be run in an embarrassingly parallel fashion.
for seed in {1..3}
do
Rscript spline_runfile.r $seed 1
Rscript fpca_runfile.r $seed 1
done
###################

# Verify convergence by plotting the traceplots
Rscript mcmc_outfile.r 1
Rscript mcmc_outfile.r 0

# This loop performs the simulation run for 100 datasets
# The for-loop above can be run in an embarrassingly parallel fashion.
for seed in {1..100}
do
Rscript spline_runfile.r $seed 0
Rscript fpca_runfile.r $seed 0
done
###################

# Run all evaluation metrics to compare the results
Rscript evaluation_metrics.r

# Plot the additional visualizations in either the paper or simulation
Rscript visualizations.r

# Plot posterior probability plots for one particular dataset
Rscript mcmc_chart_plot.r 1 1
Rscript mcmc_chart_plot.r 1 0
