source("final_routine.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

load('Data/data_format.rda')
load('Data/true_pars.rda')
load('Data/par_index.rda')

init_par = pars

init_state = data_format[,"true_state"]

prior_mean = rep(0, 3)
prior_sd = rep(20, 3)
prior_par = data.frame( prior_mean= prior_mean,
                        prior_sd= prior_sd)

temp_data = as.matrix(data_format); rownames(temp_data) = NULL
id = temp_data[,"id"]
y = temp_data[,"y"]
t = temp_data[,"t"]
steps = 10000
burnin = 5000

# Creating the B-spline.
# We will change K in a dense grid around 10 to determine a metric to test the 
# optimal number of basis splines
# Since each subject has the same observation points, each B1(i) = B2(i) = B
# Otherwise we would have to create a list for each of the different subjectss
K = 10
unique_t = unique(t) 
basis_init = create.bspline.basis(range(unique_t),norder=4,nbasis=K)
basis_init_eval_pts = getbasismatrix(unique_t, basis_init, nderiv=0)
B_1 = B_2 = basis_init_eval_pts


s_time = Sys.time()

mcmc_out = mcmc_routine(y, t, id, init_par, prior_par, par_index,
                        steps, burnin, ind, init_state, B1, B2)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, ".rda"))