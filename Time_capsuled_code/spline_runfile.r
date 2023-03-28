source("spline_routine.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

trial_num = 1

load('Data/data_format.rda')
load('Data/true_pars.rda')
load('Data/par_index.rda')
load('Data/big_B.rda')

init_par = pars

init_state = data_format[,"true_state"]

prior_mean = rep(0, 3)
prior_sd = rep(20, 3)
prior_par = data.frame( prior_mean= prior_mean,
                        prior_sd= prior_sd)

# These are initialized such that all Z_1 = Z_2 = 1
B_1_master = B_2_master = big_B[['10']]

temp_data = as.matrix(data_format); rownames(temp_data) = NULL
id = temp_data[,"id"]
y = temp_data[,"y"]
t = temp_data[,"t"]
steps = 10000
burnin = 5000


s_time = Sys.time()

mcmc_out = mcmc_routine(y, t, id, init_par, prior_par, par_index,
                        steps, burnin, ind, init_state, B_1_master, B_2_master)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, "_bspline.rda"))