# -------------------------------------------------- #
# The code below was written with the assistance of  #
# ------------- Ramsay et al. (2009) --------------- #
# -------------- Sousa et al (2023) ---------------- #
# -------------------------------------------------- #
library(fda)

set.seed(102)

data_format = NULL

# Setting up the true parameters / functional curves ---------------------------

# The same evaluation points will be used for all patients
t = seq(0.01,1,by = 0.01) # NEXT: try t_i = runif(100, 0, 1)
sigma2 = 0.1^2

# Generate a B-spline of degree 3 (i.e. order 4) with K=10 basis functions for each state
basis_s1 = create.bspline.basis(range(t),norder=4,nbasis=10)
basis_s2 = create.bspline.basis(range(t),norder=4,nbasis=10)

basis_s1_eval_pts = getbasismatrix(t, basis_s1, nderiv=0)
basis_s2_eval_pts = getbasismatrix(t, basis_s2, nderiv=0)

# Whichever beta we set to 0 is an indication that it is not necessary
coef_s1 = c(-2, 0,  1.5, 1.5, 0, -1, -0.5, -1, 0, 0) # beta_1
coef_s2 = c( 1, 3, -0.5,  -1, 0,  2,  0.5,  1, 2, 1) # beta_2

true_fnc_1 = basis_s1_eval_pts %*% coef_s1
true_fnc_2 = basis_s2_eval_pts %*% coef_s2
# plot(true_fnc_1, ylim = c(min(c(true_fnc_1, true_fnc_2)), max(c(true_fnc_1, true_fnc_2))))
# points(true_fnc_2)

pars = c( 0,       # init
          -2, -2,  # omega
          sigma2,  # sigma2
          coef_s1, # beta_1
          coef_s2) # beta_2

par_index = list( init=1, omega = 2:3, sigma2 = 4,
                  beta_1 = 5:14,
                  beta_2 = 15:24)
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

save(data_format, file = 'Data/data_format.rda')
save(pars, file = 'Data/true_pars.rda')
save(par_index, file = 'Data/par_index.rda')
save(y_mat, file = 'Data/y_mat.rda')


# Creating the B-spline to be used in the runfile
# We will change K in a dense grid around 10 to determine a metric to test the 
# optimal number of basis splines
# Since each subject has the same observation points, each B1(i) = B2(i) = B
# Otherwise we would have to create a list for each of the different subjectss
K = 5:20
big_B = vector(mode = 'list', length = length(K))

for(k in K) {
    basis_init = create.bspline.basis(range(t),norder=4,nbasis=k)
    basis_init_eval_pts = getbasismatrix(t, basis_init, nderiv=0)
    big_B[[k-4]] = basis_init_eval_pts
}
names(big_B) = as.character(K)
save(big_B, file = 'Data/big_B.rda')



