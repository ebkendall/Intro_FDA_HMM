source("final_routine.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

trial_num = 1

load('Data/data_format.rda')
load('Data/par_index.rda')
load('Data/y_mat.rda')

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

# Establishing the basis from FPCA --------------------------------------------
# B_1 = B_2 = big_B[['10']]
N = nrow(y_mat)
w = diff(t)[1]

V_mat = (1/N) * (t(y_mat) %*% y_mat)
eigen_V = eigen(V_mat)
rho = eigen_V$values * w
xi = eigen_V$vectors * (1/sqrt(w))

# percentage variability explained needs to be more than 0.02
# p_comp_keep = which(rho/sum(rho) > 0.02) 
p_comp_keep = 1:10

B_1 = B_2 = xi[,p_comp_keep]

# Tidemann miller then multiplies each eigenfunction by its respective eigenvalue

init_par = c( 0,      # init
              -2, -2,  # omega
              0.1^2,  # sigma2
              rep(0, length(p_comp_keep)), # beta_1
              rep(0, length(p_comp_keep))) # beta_2

s_time = Sys.time()

mcmc_out = mcmc_routine(y, t, id, init_par, prior_par, par_index,
                        steps, burnin, ind, init_state, B_1, B_2)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, "_fpca.rda"))