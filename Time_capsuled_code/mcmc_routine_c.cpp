#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

//  FUNCTIONS: ---------------------------------------------------------------

// [[Rcpp::export]]
arma::vec update_theta_j(const arma::vec &pars, const arma::field<arma::uvec> &par_index, const int j) {
    
    // par_index: (0) init, (1) omega, (2) sigma2, (3) beta_1, (4) beta_2
    //            (5) theta_1, (6) theta_2, (7) Z_1, (8) Z_2
    
    arma::vec theta_new(par_index(5).n_elem);
    arma::vec zed_1 = pars.elem(par_index(7) - 1);
    arma::vec zed_2 = pars.elem(par_index(8) - 1);
    
    // prior specifications for theta
    double mu_1 = 0.5;
    double mu_2 = 0.5;
    
    if(j == 1) {
        for(int i = 0; i < theta_new.n_elem; i++) {
            double alpha = mu_1 + zed_1(i);
            double beta  = 2 - mu_1 - zed_1(i);
            
            theta_new(i) = r_4beta(alpha, beta, 0, 1);
        }
    }
    
    if(j == 2) {
        for(int i = 0; i < theta_new.n_elem; i++) {
            double alpha = mu_2 + zed_2(i);
            double beta  = 2 - mu_2 - zed_2(i);
            
            theta_new(i) = r_4beta(alpha, beta, 0, 1);
        }
    }
    
    return theta_new;
}

arma::mat update_B_mat(arma::vec zed_j, arma::mat &B_j) {
    
    // par_index: (0) init, (1) omega, (2) sigma2, (3) beta_1, (4) beta_2
    //            (5) theta_1, (6) theta_2, (7) Z_1, (8) Z_2
    
    arma::mat new_B_j = B_j;
    
    for(int i = 0; i < zed_j.n_elem; i++) {
        new_B_j.col(i) = zed_j(i) * new_B_j.col(i);
    }
    
    return new_B_j;
}

// [[Rcpp::export]]
arma::field<arma::mat> update_both_B(const arma::vec &pars, const arma::field<arma::uvec> &par_index,
                                     const arma::mat &B_1, const arma::mat &B_2) {
    // par_index: (0) init, (1) omega, (2) sigma2, (3) beta_1, (4) beta_2
    //            (5) theta_1, (6) theta_2, (7) Z_1, (8) Z_2
    
    arma::vec zed_1 = pars.elem(par_index(7) - 1);
    arma::vec zed_2 = pars.elem(par_index(8) - 1);
    
    arma::mat new_B_1 = B_1;
    arma::mat new_B_2 = B_2;
    
    for(int i = 0; i < zed_1.n_elem; i++) {
        new_B_1.col(i) = zed_1(i) * new_B_1.col(i);
        new_B_2.col(i) = zed_2(i) * new_B_2.col(i);
    }
    
    arma::field<arma::mat> B_list(2);
    B_list(0) = new_B_1; 
    B_list(1) = new_B_2;
    
    return B_list;
}

// [[Rcpp::export]]
arma::vec update_beta_j(const arma::vec &pars, const arma::field<arma::uvec> &par_index,
                        const arma::field<arma::vec> &W_1, const arma::field<arma::vec> &W_2,
                        const arma::vec &y, const arma::vec &id,
                        const int j, const arma::vec &EIDs,
                        const arma::mat &B_1, const arma::mat &B_2) {
    // par_index: (0) init, (1) omega, (2) sigma2, (3) beta_1, (4) beta_2
    //            (5) theta_1, (6) theta_2, (7) Z_1, (8) Z_2
    
    arma::vec beta_1 = pars.elem(par_index(3) - 1);
    arma::vec beta_2 = pars.elem(par_index(4) - 1);
    double sigma2 = arma::as_scalar(pars.elem(par_index(2)- 1));
    double inv_sigma2 = 1/sigma2;
    
    // Defining the uninformed priors on beta_1 and beta_2
    arma::mat alpha_1(beta_1.n_elem, 1, arma::fill::zeros);
    arma::mat alpha_2(beta_2.n_elem, 1, arma::fill::zeros);
    double tau2 = 20;
    double inv_tau2 = 1/tau2;
    
    // Calculating the conjugate mean and variance
    arma::mat W_inv_sub(beta_1.n_elem, beta_1.n_elem, arma::fill::zeros);
    arma::vec V_sub(beta_1.n_elem, arma::fill::zeros);
    
    for (int ii = 0; ii < EIDs.n_elem; ii++) { 
        
        int i = EIDs(ii);
        
        arma::uvec sub_ind = arma::find(id == i);
        arma::vec y_i = y.elem(sub_ind);
        
        // Dealing with beta_1
        if (j == 1) {
            arma::mat W_mat = arma::diagmat(W_1(ii));
            arma::mat interm = B_1.t() * W_mat.t() * W_mat * B_1;
            W_inv_sub = W_inv_sub + interm;

            arma::mat not_W_mat = arma::diagmat(W_2(ii));
            arma::vec interm_v = B_1.t() * W_mat.t() * (y_i - not_W_mat * B_2 * beta_2);
            V_sub = V_sub + interm_v;
        } 
        // Dealing with beta_2
        else{
            arma::mat W_mat = arma::diagmat(W_2(ii));
            arma::mat interm = B_2.t() * W_mat.t() * W_mat * B_2;
            W_inv_sub = W_inv_sub + interm;
            
            arma::mat not_W_mat = arma::diagmat(W_1(ii));
            arma::vec interm_v = B_2.t() * W_mat.t() * (y_i - not_W_mat * B_1 * beta_1);
            V_sub = V_sub + interm_v;
        }
    }
    arma::mat W_inv = inv_tau2 * arma::eye(beta_1.n_elem, beta_1.n_elem) + inv_sigma2 * W_inv_sub;
    arma::mat W = arma::inv_sympd(W_inv);
    
    arma::mat V;
    if(j == 1) {
        V = inv_tau2 * alpha_1 + inv_sigma2 * V_sub;
    } else{
        V = inv_tau2 * alpha_2 + inv_sigma2 * V_sub;
    }
    
    return arma::mvnrnd(W * V, W, 1);
    
}

// [[Rcpp::export]]
double update_sigma2(const arma::vec &pars, const arma::field<arma::uvec> &par_index, const arma::vec &id,
                     const arma::vec &y, const arma::vec &EIDs, const arma::field<arma::vec> &W_1, 
                     const arma::field<arma::vec> &W_2, const arma::mat &B_1, const arma::mat &B_2) {
    // par_index: (0) init, (1) omega, (2) sigma2, (3) beta_1, (4) beta_2
    //            (5) theta_1, (6) theta_2, (7) Z_1, (8) Z_2
    int lambda1 = 1;
    int lambda2 = 1;
    
    int N_total = y.n_elem;
    
    arma::vec total = {0};
    
    for (int ii = 0; ii < EIDs.n_elem; ii++) { 
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(id == i);
        arma::vec y_i = y.elem(sub_ind);
        
        arma::mat W_1_i = arma::diagmat(W_1(ii));
        arma::mat W_2_i = arma::diagmat(W_2(ii));
        
        arma::vec beta_1 = pars.elem(par_index(3) - 1);
        arma::vec beta_2 = pars.elem(par_index(4) - 1);
        
        arma::mat diff = y_i - W_1_i*B_1*beta_1 - W_2_i*B_2*beta_2;
        
        total = total + diff.t() * diff;
    }
    
    double lambda1_new = lambda1 + 0.5 * N_total;
    double lambda2_new = lambda2 + 0.5 * arma::as_scalar(total);
    
    // Since C++ does not have an inv-gamma generator, we exploit the 
    // relationship between inv-gamma and inv-wishart
    double df = 2 * lambda1_new;
    arma::mat scale_mat(1,1,arma::fill::ones);
    scale_mat(0,0) = 2*lambda2_new;
    
    arma::mat final = riwish(df, scale_mat);
    
    return(arma::as_scalar(final));
}

// [[Rcpp::export]]
double fn_log_post_continuous(const arma::vec &pars, const arma::field<arma::vec> &prior_par, 
                              const arma::field<arma::uvec> &par_index, const arma::vec &y, 
                              const arma::vec &id, const arma::vec &EIDs,
                              const arma::mat &B_1, const arma::mat &B_2) {
    
    // par_index: (0) init, (1) omega, (2) sigma2, (3) beta_1, (4) beta_2
    //            (5) theta_1, (6) theta_2, (7) Z_1, (8) Z_2
    // "i" is the numeric ID number
    // "ii" is the index of the ID
    
    // Storage of the likelihood expression
    arma::vec in_vals(EIDs.n_elem, arma::fill::zeros);
    
    // Initial state probabilities
    arma::vec init_logit = {1, exp(arma::as_scalar(pars.elem(par_index(0) - 1)))};
    arma::vec init = init_logit / arma::accu(init_logit);
    
    // Transition probability matrix
    arma::vec omega = pars.elem(par_index(1) - 1);
    arma::mat P = { {1, exp(omega(0))},
                    {exp(omega(1)), 1} };
    arma::vec P_row_sums = arma::sum(P, 1);
    P = P.each_col() / P_row_sums;
    
    // Defining key terms
    arma::vec beta_1 = pars.elem(par_index(3) - 1);
    arma::vec beta_2 = pars.elem(par_index(4) - 1);
    double sigma2 = arma::as_scalar(pars.elem(par_index(2) - 1));
    double sigma = sqrt(sigma2);
    
    // Parallelized computation of the log-likelihood
    omp_set_num_threads(6);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        
        // Subsetting the data
        arma::uvec sub_ind = arma::find(id == i);
        arma::vec y_i = y.elem(sub_ind);
        
        arma::vec B_1_row = B_1.row(0).t();
        arma::vec B_2_row = B_2.row(0).t();
        double mean_1 = arma::dot(B_1_row, beta_1);
        double mean_2 = arma::dot(B_2_row, beta_2);
        
        double d_1 = arma::normpdf(y_i(0), mean_1, sigma);
        double d_2 = arma::normpdf(y_i(0), mean_2, sigma);
        
        arma::vec d_fill = {d_1, d_2};
        arma::mat D_i = arma::diagmat(d_fill);
        
        arma::mat init_transpose = init.t();
        
        arma::mat f_i = init_transpose * D_i;
        
        for(int k = 1; k < y_i.n_elem; k++) {
            
            arma::vec B_1_row = B_1.row(k).t();
            arma::vec B_2_row = B_2.row(k).t();
            double mean_1 = arma::dot(B_1_row, beta_1);
            double mean_2 = arma::dot(B_2_row, beta_2);
            
            double d_1 = arma::normpdf(y_i(k), mean_1, sigma);
            double d_2 = arma::normpdf(y_i(k), mean_2, sigma);
            
            arma::vec d_fill = {d_1, d_2};
            arma::mat D_i = arma::diagmat(d_fill);
            
            arma::mat temp = f_i * P * D_i;
            f_i = temp;
        }
        
        in_vals(ii) = log(arma::accu(f_i));
    }
    
    double in_value = arma::accu(in_vals);
    
    
    // Likelihood components from the Metropolis priors
    arma::vec p_mean = prior_par(0);
    arma::mat p_sd = arma::diagmat(prior_par(1));
    
    arma::vec init_val = pars.elem(par_index(0) - 1);
    arma::mat p_x = arma::join_cols(init_val, omega);
    
    double log_prior_dens = arma::as_scalar(dmvnorm(p_x.t(), p_mean, p_sd, true));
    
    in_value = in_value + log_prior_dens;
    
    return in_value;
}

double log_f_i_cpp(const int i, const int ii, const arma::vec &pars, const arma::field<arma::uvec> &par_index,
                   const arma::vec &y, const arma::vec &id, const arma::vec &B, 
                   const int &state_ind, const arma::mat &B_1, const arma::mat &B_2) {
    
    // par_index: (0) init, (1) omega, (2) sigma2, (3) beta_1, (4) beta_2
    //            (5) theta_1, (6) theta_2, (7) Z_1, (8) Z_2
    // "i" is the numeric ID number
    // "ii" is the index of the ID
    double in_value = 0;
    
    // Subsetting the data
    arma::uvec sub_ind = arma::find(id == i);
    arma::vec y_i = y.elem(sub_ind);
    
    // Current State sequence
    arma::mat b_i = B;
    
    // Initial state probabilities
    arma::vec init_logit = {1, exp(arma::as_scalar(pars.elem(par_index(0) - 1)))};
    arma::vec init = init_logit / arma::accu(init_logit);
    
    // Transition probability matrix
    arma::vec omega = pars.elem(par_index(1) - 1);
    arma::mat P = { {1, exp(omega(0))},
                    {exp(omega(1)), 1} };
    arma::vec P_row_sums = arma::sum(P, 1);
    P = P.each_col() / P_row_sums;
    
    // Defining key terms
    arma::vec beta_1 = pars.elem(par_index(3) - 1);
    arma::vec beta_2 = pars.elem(par_index(4) - 1);
    double sigma2 = arma::as_scalar(pars.elem(par_index(2) - 1));
    double sigma = sqrt(sigma2);
    
    // Defining the changed state
    int b_k = b_i(state_ind);
    double y_mean;
    if (b_k == 1){
        arma::vec B_1_row = B_1.row(state_ind).t();
        y_mean = arma::dot(B_1_row, beta_1);
    } 
    if (b_k == 2) {
        arma::vec B_2_row = B_2.row(state_ind).t();
        y_mean = arma::dot(B_2_row, beta_2);
    }
    
    // Likelihood evaluated at the state change point
    if(state_ind==0){
        // Initial State
        int b_k_2 = b_i(state_ind+1);
        in_value = in_value + log(init(b_k - 1)) + log(P( b_k - 1, b_k_2 - 1));
    } else if(state_ind == y_i.n_elem - 1) {
        // Final state
        int b_k_1 = b_i(state_ind-1);
        in_value = in_value + log(P( b_k_1 - 1, b_k - 1));
    } else{
        int b_k_1 = b_i(state_ind-1);
        int b_k_2 = b_i(state_ind+1);
        in_value = in_value + log(P( b_k_1 - 1, b_k - 1)) + log(P( b_k - 1, b_k_2 - 1));
    }
    
    in_value = in_value + arma::log_normpdf(y_i(state_ind), y_mean, sigma);
    
    return in_value;
}

// [[Rcpp::export]]
arma::field<arma::vec> update_b_i_cpp(const arma::vec &pars, const arma::field<arma::uvec> &par_index,
                                      const arma::vec &y, const arma::vec &id, const arma::vec &EIDs,
                                      arma::field<arma::vec> &B, const arma::mat &B_1, const arma::mat &B_2) {
    
    // par_index: (0) init, (1) omega, (2) sigma2, (3) beta_1, (4) beta_2
    //            (5) theta_1, (6) theta_2, (7) Z_1, (8) Z_2
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    arma::field<arma::vec> B_return(EIDs.n_elem);
    
    omp_set_num_threads(6) ;
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        
        // Subsetting the data
        arma::uvec sub_ind = arma::find(id == i);
        arma::vec y_i = y.elem(sub_ind);
        
        // Subsetting fields
        arma::vec B_temp = B(ii);
        
        // Update state sequence 1 at a time
        for (int k = 0; k < y_i.n_elem; k++) {
            
            arma::vec pr_B = B_temp;
            
            // Sample and update the next state
            int sampled_index = arma::randi(arma::distr_param(1, 2));
            
            pr_B(k) = sampled_index;
            
            double log_target_prev = log_f_i_cpp(i, ii, pars, par_index, y, id, B_temp, k, B_1, B_2);
            double log_target      = log_f_i_cpp(i, ii, pars, par_index, y, id,   pr_B, k, B_1, B_2);
            
            // Note that the proposal probs cancel in the MH ratio
            double diff_check = log_target - log_target_prev;
            double min_log = log(arma::randu(arma::distr_param(0,1)));
            if(diff_check > min_log){
                B_temp = pr_B;
            }
        }
        B_return(ii) = B_temp;
    }
    
    return B_return;
}


// [[Rcpp::export]] 
arma::vec update_Z_j(const arma::vec &pars, const arma::field<arma::vec> &prior_par, 
                     const arma::field<arma::uvec> &par_index, const arma::vec &y, 
                     const arma::vec &id, const arma::vec &EIDs,
                     arma::mat &B_1, arma::mat &B_2, const int j) {
    
    // par_index: (0) init, (1) omega, (2) sigma2, (3) beta_1, (4) beta_2
    //            (5) theta_1, (6) theta_2, (7) Z_1, (8) Z_2
    
    arma::vec zed_1_curr = pars.elem(par_index(7) - 1);
    arma::mat B_1_interm = update_B_mat(zed_1_curr, B_1);
    
    arma::vec zed_2_curr = pars.elem(par_index(8) - 1);
    arma::mat B_2_interm = update_B_mat(zed_2_curr, B_2);
    
    arma::vec theta_1    = pars.elem(par_index(5) - 1);
    arma::vec theta_2    = pars.elem(par_index(6) - 1);
    
    arma::vec zed_update;
    if(j == 1) zed_update = zed_1_curr;
    if(j == 2) zed_update = zed_2_curr;
    
    for(int i = 0; i < zed_update.n_elem; i++) {
        if(j == 1) {
            arma::vec z_k_0 = zed_1_curr; z_k_0(i) = 0;
            arma::vec z_k_1 = zed_1_curr; z_k_1(i) = 1;
            
            arma::mat B_k_0 = update_B_mat(z_k_0, B_1);
            arma::mat B_k_1 = update_B_mat(z_k_1, B_1);
            
            double f_0 = fn_log_post_continuous(pars, prior_par, par_index, y, id, EIDs, B_k_0, B_2_interm);
            double f_1 = fn_log_post_continuous(pars, prior_par, par_index, y, id, EIDs, B_k_1, B_2_interm);
            
            double interm = exp(f_0 - f_1);
            double p_success = theta_1(i) / ((1 - theta_1(i)) * interm + theta_1(i));
            
            double new_Z = R::rbinom(1, p_success);
            zed_update(i) = new_Z;
            zed_1_curr(i) = new_Z;
        }
        
        else if(j == 2) {
            arma::vec z_k_0 = zed_2_curr; z_k_0(i) = 0;
            arma::vec z_k_1 = zed_2_curr; z_k_1(i) = 1;
            
            arma::mat B_k_0 = update_B_mat(z_k_0, B_2);
            arma::mat B_k_1 = update_B_mat(z_k_1, B_2);
            
            double f_0 = fn_log_post_continuous(pars, prior_par, par_index, y, id, EIDs, B_1_interm, B_k_0);
            double f_1 = fn_log_post_continuous(pars, prior_par, par_index, y, id, EIDs, B_1_interm, B_k_1);
            
            double interm = exp(f_0 - f_1);
            double p_success = theta_2(i) / ((1 - theta_2(i)) * interm + theta_2(i));
            
            double new_Z = R::rbinom(1, p_success);
            zed_update(i) = new_Z;
            zed_2_curr(i) = new_Z;
        }
    }
    
    return zed_update;
}

// [[Rcpp::export]]
void test_functions(const arma::vec &pars, const arma::field<arma::uvec> &par_index) {
    
    for(int i = 0; i < 20; i++) {
        Rcpp::Rcout << arma::randi(arma::distr_param(1, 2)) << std::endl;
    }
    
}