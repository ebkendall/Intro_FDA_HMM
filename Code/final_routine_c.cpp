#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// Defining the Omega_List as a global variable when pre-compiling ----------
const arma::mat adj_mat = { {1, 1},
                            {1, 1} };


arma::field<arma::field<arma::mat>> Omega_set(const arma::mat &G) {
    int N = G.n_cols; // dimension of adj matrix
    
    arma::field<arma::mat> c(N);
    arma::field<arma::mat> b(N, N);
    arma::field<arma::mat> a(N);
    
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            b(i, j) = arma::mat(1, 2, arma::fill::zeros);
        }
    }
    
    // a -------------------------------------------------------
    for(int i = 0; i < N; i++) {
        arma::mat a_i(1, 2, arma::fill::zeros);
        arma::uvec sec_elem = arma::find(G.row(i) == 1);
        
        for(int j = 0; j < sec_elem.n_elem; j++) {
            int sec_ind = sec_elem(j);
            arma::uvec third_elem = arma::find(G.row(sec_ind) == 1);
            
            for(int k = 0; k < third_elem.n_elem; k++) {
                int third_ind = third_elem(k);
                arma::mat temp(1,2);
                temp(0,0) = sec_ind; temp(0,1) = third_ind;
                a_i = arma::join_vert(a_i, temp+1);
            }
        }
        
        a_i = a_i.rows(1, a_i.n_rows-1);
        a(i) = a_i;
    }
    
    // b -------------------------------------------------------
    for(int i = 0; i < N; i++) {
        arma::uvec sec_elem = arma::find(G.row(i) == 1);
        
        for(int j = 0; j < sec_elem.n_elem; j++) {
            int sec_ind = sec_elem(j);
            arma::uvec third_elem = arma::find(G.row(sec_ind) == 1);
            
            for(int k = 0; k < third_elem.n_elem; k++) {
                int third_ind = third_elem(k);
                arma::uvec fourth_elem = arma::find(G.row(third_ind) == 1);
                
                for(int l = 0; l < fourth_elem.n_elem; l++) {
                    int fourth_ind = fourth_elem(l);
                    arma::mat temp(1,2);
                    temp(0,0) = sec_ind; temp(0,1) = third_ind;
                    b(i, fourth_ind) = arma::join_vert(b(i, fourth_ind), temp+1);
                }
            }
        }
    }
    
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            if (b(i,j).n_rows > 1) {
                b(i, j) = b(i, j).rows(1, b(i,j).n_rows - 1);
            } else{
                arma::mat temp(1,2); temp(0,0) = -1; temp(0,1) = -1;
                b(i,j) = temp;
            }
        }
    }
    
    // c -------------------------------------------------------
    for(int i = 0; i < N; i++) {
        arma::mat c_i(1, 2, arma::fill::zeros);
        arma::uvec sec_elem = arma::find(G.col(i) == 1);
        
        for(int j = 0; j < sec_elem.n_elem; j++) {
            int sec_ind = sec_elem(j);
            arma::uvec third_elem = arma::find(G.col(sec_ind) == 1);
            
            for(int k = 0; k < third_elem.n_elem; k++) {
                int third_ind = third_elem(k);
                arma::mat temp(1,2);
                temp(0,0) = third_ind; temp(0,1) = sec_ind;
                c_i = arma::join_vert(c_i, temp+1);
            }
        }
        
        c_i = c_i.rows(1, c_i.n_rows-1);
        c(i) = c_i;
    }
    
    arma::field<arma::field<arma::mat>> Omega_List(3);
    Omega_List(0) = c; Omega_List(1) = b; Omega_List(2) = a;
    
    return Omega_List;
}

const arma::field<arma::field<arma::mat>> Omega_List_GLOBAL = Omega_set(adj_mat);

const double pi = 3.14159265358979323846;

//  FUNCTIONS: ---------------------------------------------------------------

arma::mat Omega_fun_cpp_new(const int k, const int n_i, const arma::vec &b_i) {
    
    arma::mat Omega_set;
    
    // b(k) is either 1 or 2 therefore subtract 1 for the index
    if (k == 1) {
        // () -> () -> 1-2
        Omega_set = Omega_List_GLOBAL(0)(b_i(2) - 1);
    } else if (k <= n_i - 2) {
        // 1-2 -> () -> () -> 1-2
        Omega_set = Omega_List_GLOBAL(1)(b_i(k - 2) - 1, b_i(k + 1) - 1);
    } else if (k == n_i - 1) {
        // 1-2 -> () -> ()
        Omega_set = Omega_List_GLOBAL(2)(b_i(n_i - 3) - 1);
    }
    
    return Omega_set;
    
}

// [[Rcpp::export]]
arma::vec update_beta_j(const arma::vec &pars, const arma::field<arma::uvec> &par_index,
                        const arma::field<arma::vec> &W_1, const arma::field<arma::vec> &W_2,
                        const arma::vec &y, const arma::vec &id,
                        const arma::field<arma::mat> &K, const int j, const arma::vec &EIDs,
                        const arma::mat &B_1, const arma::mat &B_2) {
    // par_index: (0) init, (1) omega, (2) sigma2, (3) beta_1, (4) beta_2
    
    arma::vec beta_1 = pars.elem(par_index(3) - 1);
    arma::vec beta_2 = pars.elem(par_index(4) - 1);
    double sigma2 = arma::as_scalar(pars.elem(par_index(2)- 1));
    double inv_sigma2 = 1/sigma2;
    
    // Defining the informed priors on beta_1 and beta_2
    // arma::mat alpha_1(beta_1.n_elem, 1, arma::fill::zeros);
    // arma::mat alpha_2(beta_2.n_elem, 1, arma::fill::zeros);
    arma::mat alpha_1 = {-2, 0,  1.5, 1.5, 0, -1, -0.5, -1, 0, 0};
    arma::mat alpha_2 = { 1, 3, -0.5,  -1, 0,  2,  0.5,  1, 2, 1};
    double tau2 = 10;
    double inv_tau2 = 1/tau2;
    
    // Calculating the conjugate mean and variance
    arma::mat W_inv_sub;
    arma::vec V_sub;
    
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
    
    arma::mat W_inv = inv_tau2 * arma::eye(B_1.n_rows, B_1.n_cols) + inv_sigma2 * W_inv_sub;
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
                     const arma::field<arma::vec> &B, const arma::vec &y, const arma::vec &EIDs) {
    // par_index: (0) init, (1) omega, (2) sigma2, (3) beta_1, (4) beta_2
    int a = 1;
    int b = 1;
    
    int N_total = y.n_elem;
    
    arma::vec total = {0};
    
    for (int ii = 0; ii < EIDs.n_elem; ii++) { 
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(id == i);
        arma::vec y_i = y.elem(sub_ind);
        arma::vec b_i = B(ii);
        
        arma::vec f1 = pars.elem(par_index(4) - 1);
        arma::vec f2 = pars.elem(par_index(5) - 1);
        
        f1.elem(arma::find(b_i == 2)).zeros();
        f2.elem(arma::find(b_i == 1)).zeros();
        
        arma::mat diff = y_i - f1 - f2;
        
        total = total + diff.t() * diff;
    }
    
    double a_new = a + 0.5 * N_total;
    double b_new = b + 0.5 * arma::as_scalar(total);
    
    // Since C++ does not have an inv-gamma generator, we exploit the 
    // relationship between inv-gamma and inv-wishart
    double df = 2 * a_new;
    arma::mat scale_mat(1,1,arma::fill::ones);
    scale_mat(0,0) = 2*b_new;
    
    arma::mat final = riwish(df, scale_mat);
    
    return(arma::as_scalar(final));
}

// [[Rcpp::export]]
double fn_log_post_continuous(const arma::vec &pars, const arma::field<arma::vec> &prior_par, 
                              const arma::field<arma::uvec> &par_index,
                              const arma::vec &y, arma::vec &x, const arma::vec &id,
                              const arma::field<arma::mat> &K, const arma::vec &EIDs) {
    
    // par_index: (0) init, (1) omega, (2) sigma2, (3) beta_1, (4) beta_2
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // Storage of the likelihood expression
    arma::vec in_vals(EIDs.n_elem, arma::fill::zeros);
    
    // Initial state probabilities
    arma::vec init_logit = {1, exp(arma::as_scalar(pars.elem(par_index(0) - 1)))};
    arma::vec init = init_logit / arma::accu(init_logit);
    
    // Transition probability matrix
    arma::vec beta = pars.elem(par_index(1) - 1);
    arma::mat P = {{1, exp(beta(0))},
    {exp(beta(1)), 1}};
    arma::vec P_row_sums = arma::sum(P, 1);
    P = P.each_col() / P_row_sums;
    
    // Defining key terms
    arma::vec f1 = pars.elem(par_index(4) - 1);
    arma::vec f2 = pars.elem(par_index(5) - 1);
    double sigma2 = arma::as_scalar(pars.elem(par_index(2) - 1));
    double sigma = sqrt(sigma2);
    
    // Parallelized computation of the log-likelihood
    // omp_set_num_threads(6);
    // # pragma omp parallel for
    
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        
        // Subsetting the data
        arma::uvec sub_ind = arma::find(id == i);
        arma::vec y_i = y.elem(sub_ind);
        
        double d_1 = arma::normpdf(y_i(0), f1(0), sigma);
        double d_2 = arma::normpdf(y_i(0), f2(0), sigma);
        
        arma::vec d_fill = {d_1, d_2};
        arma::mat D_i = arma::diagmat(d_fill);
        
        arma::mat init_transpose = init.t();
        
        arma::mat f_i = init_transpose * D_i;
        
        for(int k = 1; k < y_i.n_elem; k++) {
            
            double d_1 = arma::normpdf(y_i(k), f1(k), sigma);
            double d_2 = arma::normpdf(y_i(k), f2(k), sigma);
            
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
    arma::vec s = pars.elem(par_index(3) - 1);
    arma::mat p_x = arma::join_cols(arma::join_cols(init_val, beta),s);
    
    double log_prior_dens = arma::as_scalar(dmvnorm(p_x.t(), p_mean, p_sd, true));
    
    // Likelihood contribution from the rest
    arma::mat f1_mat = pars.elem(par_index(4) - 1);
    arma::mat f2_mat = pars.elem(par_index(5) - 1);
    
    arma::vec f_mean(f1.n_elem, arma::fill::zeros);
    arma::mat K_1 = K(0);
    arma::mat K_2 = K(1);
    
    // Optimized procedure for evaluating the likelihood -----------------------
    // double log_det_1 = arma::log_det_sympd(K_1);
    // double log_det_2 = arma::log_det_sympd(K_2);
    // arma::mat precision1 = arma::inv_sympd(K_1);
    // arma::mat precision2 = arma::inv_sympd(K_2);
    // double interm1 = arma::as_scalar(f1_mat.t() * precision1 * f1_mat);
    // double interm2 = arma::as_scalar(f2_mat.t() * precision2 * f2_mat);
    // double log_f1_dens = (f1.n_elem / -2) * log(2*pi) + -0.5 * (log_det_1 + interm1);
    // double log_f2_dens = (f1.n_elem / -2) * log(2*pi) + -0.5 * (log_det_2 + interm2);
    // -------------------------------------------------------------------------
    double log_f1_dens = arma::as_scalar(dmvnorm(f1_mat.t(), f_mean, K_1, true));
    double log_f2_dens = arma::as_scalar(dmvnorm(f2_mat.t(), f_mean, K_2, true));
    
    in_value = in_value + log_prior_dens + log_f1_dens + log_f2_dens;
    
    return in_value;
}

double log_f_i_cpp(const int i, const int ii, const arma::vec &pars, const arma::field<arma::uvec> &par_index,
                   const arma::vec &y, const arma::vec &id, const arma::vec &B) {
    
    // par_index: (0) init, (1) omega, (2) sigma2, (3) beta_1, (4) beta_2
    // "i" is the numeric EID number
    // "ii" is the index of the EID
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
    arma::vec beta = pars.elem(par_index(1) - 1);
    arma::mat P = {{1, exp(beta(0))},
    {exp(beta(1)), 1}};
    arma::vec P_row_sums = arma::sum(P, 1);
    P = P.each_col() / P_row_sums;
    
    // Defining key terms
    arma::vec f1 = pars.elem(par_index(4) - 1);
    arma::vec f2 = pars.elem(par_index(5) - 1);
    double sigma2 = arma::as_scalar(pars.elem(par_index(2) - 1));
    double sigma = sqrt(sigma2);
    
    
    // Full likelihood evaluation
    for(int w=0; w < y_i.n_elem;++w){
        if(w==0){
            int p_int = b_i(0);
            double y_mean;
            if (p_int == 1) y_mean = f1(0);
            if (p_int == 2) y_mean = f2(0);
            
            in_value = in_value + log(init(p_int - 1)) + arma::log_normpdf(y_i(0), y_mean, sigma);
        }
        else{
            int b_k_1 = b_i(w-1);
            int b_k = b_i(w);
            
            double y_mean;
            if (b_k == 1) y_mean = f1(w);
            if (b_k == 2) y_mean = f2(w);
            
            in_value = in_value + log(P( b_k_1 - 1, b_k - 1)) + arma::log_normpdf(y_i(w), y_mean, sigma);
        }
    }
    
    return in_value;
}

// [[Rcpp::export]]
arma::field<arma::vec> update_b_i_cpp(const arma::vec &pars, const arma::field<arma::uvec> &par_index,
                                      const arma::vec &y, const arma::vec &id, const arma::vec &EIDs,
                                      arma::field<arma::vec> &B) {
    
    // par_index: (0) init, (1) omega, (2) sigma2, (3) beta_1, (4) beta_2
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    arma::field<arma::vec> B_return(EIDs.n_elem);
    
    // omp_set_num_threads(8) ;
    // # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        
        // Subsetting the data
        arma::uvec sub_ind = arma::find(id == i);
        arma::vec y_i = y.elem(sub_ind);
        
        // Subsetting fields
        arma::vec B_temp = B(ii);
        
        // The first state is always S1, therefore we start at 1 instead of 0
        for (int k = 0; k < y_i.n_elem - 1; k++) {
            
            // arma::vec t_pts;
            // if (k == n_i - 2) {
            //     t_pts = arma::linspace(k+1, k+2, 2);
            // } else {
            //     t_pts = arma::linspace(k+1, k+3, 3);
            // }
            arma::vec t_pts = {k, k+1};
            arma::vec pr_B = B_temp;
            
            // Sample and update the two neighboring states
            arma::mat Omega_set = Omega_fun_cpp_new(k + 1, y_i.n_elem, B_temp);
            
            int sampled_index = arma::randi(arma::distr_param(1, Omega_set.n_rows));
            
            pr_B.rows(k, k+1) = Omega_set.row(sampled_index-1).t();
            
            double log_target_prev = log_f_i_cpp(i, ii, pars, par_index, y, id, B_temp);
            double log_target      = log_f_i_cpp(i, ii, pars, par_index, y, id, pr_B);
            
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
void test_functions(const arma::vec &pars, const arma::field<arma::uvec> &par_index) {
    
    int N = 2;
    Rcpp::Rcout << "Case (c) Full" << std::endl;
    for(int w=0; w < N; w++) {
        Rcpp::Rcout << "() -> () -> " << w+1 << std::endl;
        Rcpp::Rcout << Omega_List_GLOBAL(0)(w) << std::endl;
    }
    
    Rcpp::Rcout << "Case (b) Full" << std::endl;
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            Rcpp::Rcout << i+1 << "-->" << j+1 << std::endl;
            Rcpp::Rcout << Omega_List_GLOBAL(1)(i, j) << std::endl;
        }
    }
    
    Rcpp::Rcout << "Case (a) Full" << std::endl;
    for(int w=0; w < N; w++) {
        Rcpp::Rcout << w + 1 << " -> () -> ()" << std::endl;
        Rcpp::Rcout << Omega_List_GLOBAL(2)(w) << std::endl;
    }
    
}