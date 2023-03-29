# -------------------------------------------------- #
# The code below was written with the assistance of  #
# ------------- Ramsay et al. (2009) --------------- #
# -------------- Sousa et al (2023) ---------------- #
# -------------------------------------------------- #
library(fda)

set.seed(2022)

# Setting up the true parameters / functional curves ---------------------------

# The same evaluation points will be used for all patients
t = seq(0.01,1,by = 0.01) # NEXT: try t_i = runif(100, 0, 1)
sigma2 = 0.1^2

# Generate a B-spline of degree 3 (i.e. order 4) with K=10 basis functions for each state
basis_s1 = create.bspline.basis(range(t),norder=4,nbasis=10)
basis_s2 = create.bspline.basis(range(t),norder=4,nbasis=10)

basis_s1_eval_pts = getbasismatrix(t, basis_s1, nderiv=0)
basis_s2_eval_pts = getbasismatrix(t, basis_s2, nderiv=0)

# Generate the beta coefficients (not random effects yet)
beta_1 = c(-2, -4.5, -1.5, -3.5, -4.5, -1, -0.5, -1, -4, -1) 
beta_2 = c( 3,    6,  3.5,    1,  2.5,  5,  0.5,  2,  3,  7) 

# Generate the true Z_1 and Z_2 (not random effects yet)
theta_1 = rep(0.75, 10)
theta_2 = rep(0.75, 10)
Z_1 = Z_2 = NULL
for(i in 1:10) {
    Z_1 = c(Z_1, rbinom(n = 1, prob = theta_1[i], size = 1))
    Z_2 = c(Z_2, rbinom(n = 1, prob = theta_2[i], size = 1))
}

true_fnc_1 = basis_s1_eval_pts %*% (beta_1 * Z_1)
true_fnc_2 = basis_s2_eval_pts %*% (beta_2 * Z_2)
# plot(true_fnc_1, ylim = c(min(c(true_fnc_1, true_fnc_2)), max(c(true_fnc_1, true_fnc_2))))
# points(true_fnc_2)
# abline(h = 0)


pars = c( 0, c(-4, -4), sigma2, beta_1, beta_2, theta_1, theta_2, Z_1, Z_2)  

par_index = list( init=1, omega = 2:3, sigma2 = 4, beta_1 = 5:14, beta_2 = 15:24,
                  theta_1 = 25:34, theta_2 = 35:44, Z_1 = 45:54, Z_2 = 55:64)

save(pars, file = 'Data/true_pars.rda')
save(par_index, file = 'Data/true_par_index.rda')

fnc_vals = list('true_fnc_1' = true_fnc_1,
                'true_fnc_2' = true_fnc_2)
save(fnc_vals, file = 'Data/true_fnc_vals.rda')


for(seed in 1:100) {
    
    print(seed)
    set.seed(seed)
    
    data_format = NULL
    y_mat = matrix(nrow = 50, ncol = length(t))
    
    for(ind in 1:50){
        
        # Generate the s_i ------------------------------------------------------------
        # s_i time homogeneous Markov process
        init_prob = c(1, exp(pars[par_index$init]))
        init_prob = init_prob / sum(init_prob)
        
        tpm = matrix(c(1, exp(pars[par_index$omega][1]),
                       exp(pars[par_index$omega][2]), 1), nrow=2, byrow = T)
        tpm = tpm / rowSums(tpm)
        
        s = sample(x = 1:2, size = 1, prob = init_prob)
        for(i in 2:length(t)){
            s = c(s, sample(x = 1:2, size = 1, prob = tpm[tail(s, 1), ]))
        }
        
        # Generate the data -----------------------------------------------------------
        y = NULL
        
        for(i in 1:length(t)) {
            if(s[i] == 1) {
                y = c(y, true_fnc_1[i] + rnorm(1, mean = 0, sd = sqrt(pars[par_index$sigma2])))
            } else {
                y = c(y, true_fnc_2[i] + rnorm(1, mean = 0, sd = sqrt(pars[par_index$sigma2])))
            }
        }
        y_mat[ind, ] = y
        temp = cbind(ind, y, t, s)
        data_format = rbind(data_format, temp)
    }
    colnames(data_format) = c('id', 'y', 't', 'true_state')
    
    save(data_format, file = paste0('Data/data_format_', seed, '.rda'))
    save(y_mat, file = paste0('Data/y_mat_', seed, '.rda'))
}

# Creating the B-splines to be used in the runfile
# We will change K in a dense grid around 10 to determine the optimal K
# Since each subject has the same # of observation points, each B1(i) = B2(i) = B
K = c(5, 8, 9, 10, 11, 12, 15, 20, 25)
big_B = vector(mode = 'list', length = length(K))

for(k in 1:length(K)) {
    basis_init = create.bspline.basis(range(t),norder=4,nbasis=K[k])
    basis_init_eval_pts = getbasismatrix(t, basis_init, nderiv=0)
    big_B[[k]] = basis_init_eval_pts
}

names(big_B) = as.character(K)
save(big_B, file = 'Data/big_B.rda')



