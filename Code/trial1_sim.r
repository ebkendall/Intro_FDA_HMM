library(mvtnorm)

set.seed(100)

data_format = NULL

# Setting up the true parameters / functional curves ---------------------------
x = seq(1, 100, by = 0.5)

s1 = 28
s2 = 38

U_1 = 1/(s1*sqrt(2 * pi))
U_2 = 1/(s2*sqrt(2 * pi))

# covariance for f_1 and f_2 
# U_j exp(-(x-t)^2  /  2s_j^2)
K_1 = K_2 = matrix(nrow=length(x), ncol = length(x))

rownum = 1
for(i in x) {
    colnum = 1
    for(j in x) {
        K_1[rownum, colnum] = U_1 * exp(-((i - j)^2) / (2 * s1^2))
        K_2[rownum, colnum] = U_2 * exp(-((i - j)^2) / (2 * s2^2))
        colnum = colnum + 1
    }
    rownum = rownum + 1
}

f_1 = rmvnorm(n = 1, mean = rep(0, length(x)), sigma = K_1)
f_2 = rmvnorm(n = 1, mean = rep(0, length(x)), sigma = K_2)
plot(x, f_1, type = 'l', ylim = c(min(min(f_1), min(f_2)), max(max(f_1), max(f_2))))
lines(x, f_2, col = 'red')

pars = c( 0,       # init
          -2, -2,  # beta
          0.00005, # sigma2
          s1, s2 , # s1 and s2
          f_1, f_2) 

# par_index = list( init=1, beta = 2:3, sigma2 = 4, s = 5:6,
#                   f1 = 7:205,
#                   f2 = 206:length(pars))
par_index = list( init=1, beta = 2:3, sigma2 = 4, s = 5:6,
                  f1 = 7:106,
                  f2 = 107:length(pars))

for(ind in 1:100){
    
    # Generate the z_i ------------------------------------------------------------
    # z_i iid
    # z_1 = sample(x = 1:2, size = length(x), prob = c(0.7,0.3), replace = T)
    
    # z_i markov
    init_prob = c(1, exp(pars[par_index$init]))
    init_prob = init_prob / sum(init_prob)
    
    tpm = matrix(c(1, exp(pars[par_index$beta][1]),
                   exp(pars[par_index$beta][2]), 1), nrow=2, byrow = T)
    tpm = tpm / rowSums(tpm)
    
    z_2 = sample(x = 1:2, size = 1, prob = init_prob)
    for(i in 2:length(x)){
        z_2 = c(z_2, sample(x = 1:2, size = 1, prob = tpm[tail(z_2, 1), ]))
    }
    
    # Generate the data -----------------------------------------------------------
    
    y_2 = NULL
    
    for(i in 1:length(x)) {
        # if(z_1[i] == 1) y_1 = c(y_1, f_1[i])
        # else y_1 = c(y_1, f_2[i])
        if(z_2[i] == 1) y_2 = c(y_2, pars[par_index$f1][i])
        else y_2 = c(y_2, pars[par_index$f2][i])
    }
    
    sigma2 = pars[par_index$sigma2]
    
    y_2 = y_2 + sqrt(sigma2) * rnorm(length(y_2))
    # points(x, y_2)
    
    temp = cbind(ind, y_2, x, z_2)
    data_format = rbind(data_format, temp)
}

colnames(data_format) = c('id', 'y', 'x', 'true_state')

save(data_format, file = 'Data/data_format.rda')
save(pars, file = 'Data/true_pars.rda')
save(par_index, file = 'Data/par_index.rda')



# Initial test of getting splines to fit data
# library(splines)
# 
# deg   = 3
# ndx   = deg * length(x)^(1/5) 
# xr    = max(x)
# xl    = min(x)
# xmax  = xr + 0.01 * (xr - xl)
# xmin  = xl - 0.01 * (xr - xl)
# dt    = (xmax - xmin) / ndx
# knots = seq(xmin - deg * dt, xmax + deg * dt, by = dt)
# B     = splineDesign(knots = knots, x = x, ord = deg + 1, derivs = 0,outer.ok = TRUE)
# 
# 
# # Fit model
# fit = lsfit(B, y_1, intercept = T)
# lines(x, cbind(1, B) %*% fit$coefficients, col = 2)



