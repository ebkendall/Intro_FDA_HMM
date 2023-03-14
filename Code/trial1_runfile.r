source("trial1_routine.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

load('Data/data_format.rda')
load('Data/true_pars.rda')
load('Data/par_index.rda')

init_par = pars

init_state = data_format[,"true_state"]

prior_mean = rep(0, 5)
prior_sd = rep(5, 5)
prior_par = data.frame( prior_mean= prior_mean,
                        prior_sd= prior_sd)

temp_data = as.matrix(data_format); rownames(temp_data) = NULL
id = temp_data[,"id"]
y = temp_data[,"y"]
x = temp_data[,"x"]
steps = 30000
burnin = 5000

s_time = Sys.time()

mcmc_out = mcmc_routine(y, x, id, init_par, prior_par, par_index,
                        steps, burnin, n_cores, ind, init_state)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, ".rda"))