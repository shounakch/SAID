#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
//Rcpp::plugins(openmp)

using namespace arma;

// [[Rcpp::export]]
double random_gamma(double a) {
  return R::rgamma(a, 1.0);
}

// [[Rcpp::export]]
Rcpp::List rejection_sampler(arma::mat X1,
                            arma::mat X2,
                            double pen_param,
                            double delta_sq,
                            double tau1sq,
                            double tau2sq,
                            arma::mat cov1,
                            arma::mat cov2,
                            arma::mat prec1,
                            arma::mat prec2){
  
  //Returns the log(1 - exp(-c(y) * kappa)) for rejected samples y.
  
  int mask = 0;
  int M = X1.n_cols;
  int rej_len = 0;
  arma::vec penalty_stor(1, fill::zeros);
  penalty_stor(0) = -1;                 //First element, to be discarded later
  arma::vec mean_vec(M, fill::zeros);
  
  arma::mat C11 = delta_sq * tau1sq * cov1;
  arma::mat C12 = delta_sq * tau1sq * cov2;
  arma::mat C21 = delta_sq * tau2sq * cov1;
  arma::mat C22 = delta_sq * tau2sq * cov2;
  
  mask = 0;                           //Define mask variable for rejection sampling
  
  double QF13 = 0;
  double QF24 = 0;
  
  while((mask == 0) & (rej_len <= 1000)){
    
    arma::vec theta1 = mvnrnd(mean_vec, C11, 1);
    arma::vec theta2 = mvnrnd(mean_vec, C21, 1);
    arma::vec phi1 = mvnrnd(mean_vec, C12, 1);
    arma::vec phi2 = mvnrnd(mean_vec, C22, 1);
    
    arma::vec Pfn1 = square(X1 * theta1);
    arma::vec Pfn2 = square(X2 * phi1);
    arma::vec Pfn = Pfn1 % Pfn2;
    
    arma::vec Nfn1 = square(X1 * theta2);
    arma::vec Nfn2 = square(X2 * phi2);
    arma::vec Nfn = Nfn1 % Nfn2;
    
    double penalty_term = mean(Pfn) * mean(Nfn);
    
    double unif_rv = randu();
    
    if(unif_rv > exp(-pen_param * penalty_term)){
      
      arma::vec vek(1, fill::zeros);
      vek(0) = penalty_term;
      
      penalty_stor = join_cols(penalty_stor, vek);
      
      double QF1 = as_scalar(theta1.t() * (prec1 * theta1));
      double QF2 = as_scalar(theta2.t() * (prec1 * theta2));
      double QF3 = as_scalar(phi1.t() * (prec2 * phi1));
      double QF4 = as_scalar(phi2.t() * (prec2 * phi2));
      
      QF13 = QF13 + (QF1 + QF3);
      QF24 = QF24 + (QF2 + QF4);
      
      rej_len = rej_len + 1;
      
      mask = 0;
      
    }else{
      
      mask = 1;
      
    }
    
    // Rcpp::Rcout << "Rejection Length = " << rej_len << std::endl;
    
  }
  
  // Rcpp::Rcout << "Rejection Fine!" << std::endl;
  
  //Length of penalty_stor = 1 if no rejections
  //Otherwise, length of penalty_stor >= 2.
  
  //return penalty_stor;      //Has -1 at the beginning
  
  return Rcpp::List::create(Rcpp::Named("penalty_stor") = penalty_stor,
                     Rcpp::Named("rej_len") = rej_len,
                     Rcpp::Named("QF13") = QF13,
                     Rcpp::Named("QF24") = QF24);
  
}

// [[Rcpp::export]]
double sigmasq_sampler(arma::vec R, int n) {
  
  double sigma_sq = 1.0;
  
  double beta_sigsq = accu(square(R))/2;
  double alpha_sigsq = n / 2;
  
  sigma_sq = beta_sigsq / random_gamma(alpha_sigsq);
  return sigma_sq;
  
}

// [[Rcpp::export]]
arma::vec maineffects_sampler(arma::vec R, 
                              arma::mat X, 
                              arma::mat Psi_inv, 
                              double sigma_sq){
  
  arma::mat M1 = ((X.t()) * X)/sigma_sq;
  arma::mat M2 = chol(Psi_inv + M1);
  arma::mat M2inv = inv(M2);
  arma::mat betaVar = M2inv * (M2inv.t());
  
  arma::vec betaMean = betaVar * ((X.t() * R)/sigma_sq);
  
  arma::mat sampled_ME_coeff = mvnrnd(betaMean, betaVar, 1);
  return sampled_ME_coeff.col(0);
  
  //Make this (much) faster using Rue's algorithm.
  
}

// [[Rcpp::export]]
double pot_MALA(arma::vec R, 
                arma::mat X1, 
                arma::mat X2, 
                arma::vec param, 
                arma::mat S1, 
                arma::mat S2, 
                double sigma_sq,
                double delta_sq,
                Rcpp::List rej_obj){
  
  int M = X1.n_cols;
  int n = X1.n_rows;
  arma::vec rej_pen = rej_obj["penalty_stor"];
  int l = rej_pen.n_elem;
  
  arma::vec theta1 = param(span(0, M-1));
  arma::vec theta2 = param(span(M, (2*M)- 1));
  arma::vec phi1 = param(span((2*M), (3*M)-1));
  arma::vec phi2 = param(span((3*M), (4*M) - 1));
  
  double kappa = param(4*M);
  double pen_param = exp(kappa);
  
  double logtau1 = param(4*M + 1);
  double tau1sq = exp(2 * logtau1);
  double logtau2 = param(4*M + 2);
  double tau2sq = exp(2 * logtau2);
  
  arma::vec mainpart11 = X1 * theta1;
  arma::vec mainpart12 = X2 * phi1;
  arma::vec mainpart21 = X1 * theta2;
  arma::vec mainpart22 = X2 * phi2;
  
  arma::vec Pfn1 = square(mainpart11);
  arma::vec Pfn2 = square(mainpart12);
  arma::vec Pfn = Pfn1 % Pfn2;
  
  arma::vec Nfn1 = square(mainpart21);
  arma::vec Nfn2 = square(mainpart22);
  arma::vec Nfn = Nfn1 % Nfn2;
  
  arma::vec h = Pfn - Nfn;
  
  double pot_lik = accu(square(R - h))/(2*sigma_sq);
  
  arma::mat c_11mat = 0.5*(theta1.t() * (S1 * theta1));
  double c_11 = c_11mat(0,0);
  arma::mat c_21mat = 0.5*(theta2.t() * (S1 * theta2));
  double c_21 = c_21mat(0,0);
  arma::mat c_12mat = 0.5*(phi1.t() * (S2 * phi1));
  double c_12 = c_12mat(0,0);
  arma::mat c_22mat = 0.5*(phi2.t() * (S2 * phi2));
  double c_22 = c_22mat(0,0);
  
  double penalty_term = mean(Pfn) * mean(Nfn);
  
  double pot_rej = 0;
  if(l >= 2)
  {
    
    arma::vec new_rej_pen = rej_pen.subvec(1, l-1);
    arma::vec one_vec = ones(l-1);
    
    pot_rej = -accu(log(one_vec - exp(-pen_param * new_rej_pen)));
    
  }
  
  double QF13 = rej_obj["QF13"];
  double QF24 = rej_obj["QF24"];
  double rej_len = rej_obj["rej_len"];
  
  double pot_prior = ((c_11+c_12) / (tau1sq * delta_sq)) +    //PE from p(theta1, phi1)
    ((c_21+c_22) / (tau2sq * delta_sq)) +                     //PE from p(theta2, phi2)
    ((2 * M * logtau1) - (logtau1 - log(1 + tau1sq))) +       //PE from p(tau1^2)
    ((2 * M * logtau2) - (logtau2 - log(1 + tau2sq))) +       //PE from p(tau2^2)
    (pen_param * penalty_term) +                              //PE from penalty
    (0.5 * pow(kappa, 2.0)) +                                 //PE from prior on penalty; log(penalty) ~ N(0,1)
    (pot_rej) +                                               //PE from rejected samples for kappa
    ((2 * M * rej_len * logtau1) - (logtau1 - log(1 + tau1sq))) +       //PE from p(tau1^2), rejected samples
    ((2 * M * rej_len * logtau2) - (logtau2 - log(1 + tau2sq))) +       //PE from p(tau2^2), rejected samples
    (0.5*QF13 / (tau1sq * delta_sq)) +                                     //PE from p(rej_1, rej_3), rejected samples
    (0.5*QF24 / (tau2sq * delta_sq));                                    //PE from p(rej_2, rej_4), rejected samples
      
  double pot_total = pot_lik + pot_prior;
  
  // Rcpp::Rcout << "PE fine" << std::endl;
  
  return pot_total;
  
}

// [[Rcpp::export]]
arma::vec grad_MALA(arma::vec R, 
                    arma::mat X1, 
                    arma::mat X2, 
                    arma::vec param,
                    arma::mat S1, 
                    arma::mat S2, 
                    double sigma_sq,
                    double delta_sq,
                    Rcpp::List rej_obj){
  
  int M = X1.n_cols;
  int n = X1.n_rows;
  arma::vec rej_pen = rej_obj["penalty_stor"];
  int l = rej_pen.n_elem;
  
  arma::vec theta1 = param(span(0, M-1));
  arma::vec theta2 = param(span(M, (2*M)- 1));
  arma::vec phi1 = param(span((2*M), (3*M)-1));
  arma::vec phi2 = param(span((3*M), (4*M) - 1));
  
  double kappa = param(4*M);
  double pen_param = exp(kappa);
  
  double logtau1 = param(4*M + 1);
  double tau1sq = exp(2 * logtau1);
  double logtau2 = param(4*M + 2);
  double tau2sq = exp(2 * logtau2);
  
  arma::vec mainpart11 = X1 * theta1;
  arma::vec mainpart12 = X2 * phi1;
  arma::vec mainpart21 = X1 * theta2;
  arma::vec mainpart22 = X2 * phi2;
  
  arma::vec Pfn1 = square(mainpart11);
  arma::vec Pfn2 = square(mainpart12);
  arma::vec Pfn = Pfn1 % Pfn2;
  
  arma::vec Nfn1 = square(mainpart21);
  arma::vec Nfn2 = square(mainpart22);
  arma::vec Nfn = Nfn1 % Nfn2;
  
  arma::vec h = Pfn - Nfn;
  
  // Gradient from likelihood
  
  arma::vec v = h - R;
  
  arma::vec grad_lik_part11 = X1.t() * (Pfn2 % (mainpart11 % v));
  arma::vec grad_lik_part21 = - X1.t() * (Nfn2 % (mainpart21 % v));
  arma::vec grad_lik_part12 = X2.t() * (Pfn1 % (mainpart12 % v));
  arma::vec grad_lik_part22 = - X1.t() * (Nfn1 % (mainpart22 % v));
  
  arma::vec grad_lik = (2/sigma_sq)*join_cols(grad_lik_part11,
                        grad_lik_part21,
                        grad_lik_part12,
                        grad_lik_part22);
  
  arma::vec grad_prior = join_cols((S1 * theta1) / (tau1sq * delta_sq),
                                   (S1 * theta2) / (tau2sq * delta_sq),
                                   (S2 * phi1) / (tau1sq * delta_sq),
                                   (S2 * phi2) / (tau2sq * delta_sq));
  
  double penalty_term = mean(Pfn) * mean(Nfn);
  
  double int_P1 = mean(Pfn1);
  double int_P2 = mean(Pfn2);
  double int_N1 = mean(Nfn1);
  double int_N2 = mean(Nfn2);
  
  arma::mat mat1 = (X1.t() * X1) / n ;
  arma::mat mat2 = (X2.t() * X2) / n ;
  
  arma::vec grad_prior_penalty = 2 * join_cols((int_P2 * int_N1 * int_N2) * (mat1 * theta1),
                                               (int_P1 * int_P2 * int_N2) * (mat1 * theta2),
                                               (int_P1 * int_N1 * int_N2) * (mat2 * phi1),
                                               (int_P1 * int_P2 * int_N1) * (mat2 * phi2));
  
  arma::vec grad_total = grad_lik + grad_prior + (pen_param * grad_prior_penalty);
  
  // Construct gradient wrt pen_param
  
  double pot_grad_rej = 0;
  if(l >= 2)
  {
    
    arma::vec new_rej_pen = rej_pen.subvec(1, l-1);
    arma::vec one_vec = ones(l-1);
    pot_grad_rej = -pen_param * accu(new_rej_pen / (exp(pen_param * new_rej_pen) - one_vec));
    
  }
  
  double grad_pen = (pen_param * penalty_term) +            //Penalty
    (log(pen_param)) +                      //Prior on penalty  
    (pot_grad_rej); //Rejection sampling
  
  arma::vec grad_pen_vec(1, fill::zeros);
  grad_pen_vec(0) = grad_pen;
  
  ///// Construct grad_tau1 and grad_tau2
  
  arma::mat c_11mat = (theta1.t() * (S1 * theta1));
  double c_11 = c_11mat(0,0);
  arma::mat c_21mat = (theta2.t() * (S1 * theta2));
  double c_21 = c_21mat(0,0);
  arma::mat c_12mat = (phi1.t() * (S2 * phi1));
  double c_12 = c_12mat(0,0);
  arma::mat c_22mat = (phi2.t() * (S2 * phi2));
  double c_22 = c_22mat(0,0);
  
  double QF13 = rej_obj["QF13"];
  double QF24 = rej_obj["QF24"];
  int rej_len = rej_obj["rej_len"];
  
  double grad_tau1 = (2 * M * (1 + rej_len)) - 
    ((c_11 + c_12 + QF13) / (tau1sq * delta_sq)) - ((1 - tau1sq) / (1 + tau1sq));
  
  arma::vec grad_tau1_vec(1, fill::zeros);
  grad_tau1_vec(0) = grad_tau1;
  
  double grad_tau2 = (2 * M * (1 + rej_len)) 
    - ((c_21 + c_22 + QF24) / (tau2sq * delta_sq)) - ((1 - tau2sq) / (1 + tau2sq));
  
  arma::vec grad_tau2_vec(1, fill::zeros);
  grad_tau2_vec(0) = grad_tau2;
  
  // Rcpp::Rcout << "Gradient Fine" << std::endl;
  
  return join_cols(grad_total, grad_pen_vec, grad_tau1_vec, grad_tau2_vec);
  
}

//[[Rcpp::export]]
Rcpp::List sq_sampler(arma::vec R, 
                      arma::mat X1, 
                      arma::mat X2, 
                      arma::mat S1, 
                      arma::mat S2, 
                      double sigma_sq,
                      double delta_sq,
                      arma::vec old_param, 
                      double eps_MALA,
                      arma::mat precond_mat,
                      arma::mat precond_mat_inv,
                      int L_HMC,
                      Rcpp::List rej_obj){
  
  int M = old_param.n_elem;
  
  arma::vec new_param = old_param;
  arma::mat sig1;
  sig1 = precond_mat;
  arma::vec new_rho = mvnrnd(vec(M, fill::zeros), sig1, 1);
  
  arma::vec current_param = new_param;
  arma::vec current_rho = new_rho;
  
  //Half step of momentum
  
  new_rho = new_rho - (0.5 * eps_MALA * grad_MALA(R, X1, X2,
                                                  new_param, 
                                                  S1, S2, sigma_sq, delta_sq,
                                                  rej_obj));
  
  for(int l_ind=0; l_ind < L_HMC; ++l_ind){
    
    //Full position step
    
    // new_param = new_param + ((eps_MALA / pow(c_HMC, 2.0)) * new_rho);
    
    new_param = new_param + (eps_MALA * (precond_mat_inv * new_rho));
    
    if(l_ind != L_HMC - 1){
      
      //Full Momentum Step
      
      new_rho = new_rho - (eps_MALA * grad_MALA(R, X1, X2,
                                                new_param,
                                                S1, S2, sigma_sq, delta_sq,
                                                rej_obj));
      
    }
    
  }
  
  //Half momentum step
  
  new_rho = new_rho - (0.5 * eps_MALA * grad_MALA(R, X1, X2,
                                                  new_param, 
                                                  S1, S2, sigma_sq, delta_sq,
                                                  rej_obj));
  
  new_rho = -new_rho; //To maintain reversibility
  
  double current_U = pot_MALA(R, X1, X2, current_param,  
                              S1, S2, sigma_sq, delta_sq, rej_obj);
  
  arma::mat current_K_mat = 0.5 * (current_rho.t() * (precond_mat_inv * current_rho));
  // double current_K = accu(square(current_rho))/(2 * pow(c_HMC, 2.0));
  double current_K = current_K_mat(0,0);
  
  double new_U = pot_MALA(R, X1, X2, new_param,  
                          S1, S2, sigma_sq, delta_sq, rej_obj);
  
  arma::mat new_K_mat = 0.5 * (new_rho.t() * (precond_mat_inv * new_rho));
  // double new_K = accu(square(new_rho))/(2 * pow(c_HMC, 2.0));
  double new_K = new_K_mat(0,0);
  
  double U1 = randu();
  double energy_diff1 = exp(current_U - new_U + current_K - new_K);
  int sample_accept = 0;
  arma::vec param_accept(M, fill::zeros);
  
  arma::vec ediff(1, fill::zeros);
  ediff(0) = energy_diff1;
  
  if(ediff.is_finite() == TRUE){
    
    if(U1 <= energy_diff1){
      
      param_accept = new_param;
      sample_accept = 1;
      
    }
    else{
      
      param_accept = old_param;
      sample_accept = 0;
      
    }
    //print("No Issues with HMC.")
  }
  else{
    
    param_accept = old_param;
    sample_accept = 0;
    
    //print("HMC Diverged.")
  }
  
  return Rcpp::List::create(Rcpp::Named("sampled_param") = param_accept,
                            Rcpp::Named("accept")      = sample_accept);
  
}

// [[Rcpp::export]]
Rcpp::List SIDsampler_draws_adaptive_optimized(arma::vec y,
                                               arma::mat ME_mat, 
                                               arma::cube IE_list,
                                               arma::vec eps_MALA, 
                                               double c_HMC, 
                                               int L_HMC, 
                                               int MC,
                                               int n, 
                                               int p, 
                                               int p_cov,
                                               arma::mat SigmaME, 
                                               arma::mat SigmaME_inv,
                                               arma::mat SigmaInt, 
                                               arma::mat SigmaInt_inv,
                                               int ME_nspl, 
                                               int IE_nspl,
                                               int cutoff,
                                               arma::mat map_k_to_uv,
                                               arma::vec zero_ind,
                                               double accept_low,
                                               double accept_high,
                                               double accept_scale,
                                               double a_lamb,
                                               double b_lamb,
                                               Rcpp::List init_values,
                                               int precond){
  //Define storage matrices
  
  arma::vec alpha_stor(MC, fill::zeros);
  arma::vec sigmasq_stor(MC, fill::ones);
  
  arma::mat ME_coeff_stor(MC, ME_nspl*p, fill::zeros);
  arma::mat ME_scale_stor(MC, p, fill::ones);
  // arma::mat ME_scale_aux(MC, p, fill::ones);
  
  int K = p*(p-1)/2;
  // int IE_nbasis = IE_nspl * IE_nspl;
  
  arma::cube IE_pos_theta1(MC, IE_nspl, K, fill::zeros);
  arma::cube IE_neg_theta2(MC, IE_nspl, K, fill::zeros);
  arma::cube IE_pos_phi1(MC, IE_nspl, K, fill::zeros);
  arma::cube IE_neg_phi2(MC, IE_nspl, K, fill::zeros);
  
  arma::mat IE_scale_tausq1(MC, K, fill::ones);
  arma::mat IE_scale_tausq2(MC, K, fill::ones);
  arma::mat IE_scale_a(MC, K, fill::ones);
  arma::mat IE_scale_b(MC, K, fill::ones);
  
  // IE_scale_tausq1.row(0) = 0.01 * IE_scale_tausq1.row(0);
  // IE_scale_tausq2.row(0) = 0.01 * IE_scale_tausq2.row(0);
  
  arma::vec IE_scale_deltasq(MC, fill::ones);
  arma::vec IE_scale_nu(MC, fill::ones);
  
  arma::mat IE_pen(MC, K, fill::ones);
  arma::mat IE_pen_stor(MC, K, fill::ones);
  
  arma::mat all_interactions(n, K, fill::zeros);
  
  // Define preconditioning matrix storage
  
  int n_HMC_param = 4*IE_nspl + 3;
  
  arma::cube precond_mat_stor(n_HMC_param, n_HMC_param, K, fill::zeros);
  arma::cube precond_mat_inv_stor(n_HMC_param, n_HMC_param, K, fill::zeros);
  arma::cube grad_sample_stor(MC, n_HMC_param, K, fill::zeros);
  
  if(p_cov > 0){
    
    arma::mat cov_effect_stor(MC, p_cov, fill::zeros);
    
    //// Initialize ////
    
    // Intercept and sigma^2
    
    alpha_stor(0) = init_values["intercept_est"];
    sigmasq_stor(0) = init_values["sigmasq_est"];
    
    // Covariate effects (remove if p_cov = 0)
    
    arma::vec cov_effect_init_vec = init_values["cov_effect_est"];
    
    cov_effect_stor.row(0) = cov_effect_init_vec.t();
    
    // Main effects
    
    arma::vec ME_coeff_init_vec = init_values["ME_coeff_est"];
    
    ME_coeff_stor.row(0) = ME_coeff_init_vec.t();
    
    // Interaction effects
    
    arma::mat pos1_init_mat = init_values["IE_pos_part1_coeff_est"];
    arma::mat pos2_init_mat = init_values["IE_pos_part2_coeff_est"];
    arma::mat neg1_init_mat = init_values["IE_neg_part1_coeff_est"];
    arma::mat neg2_init_mat = init_values["IE_neg_part2_coeff_est"];
    
    IE_scale_deltasq(0) = 0.01;
    
    for(int k=0; k<K; ++k){
      
      // Obtain inverse index (u,v)
      
      int u = map_k_to_uv(k,1);
      int v = map_k_to_uv(k,2);
      
      (IE_pos_theta1.slice(k)).row(0) = pos1_init_mat.col(k).t();
      (IE_neg_theta2.slice(k)).row(0) = neg1_init_mat.col(k).t();
      (IE_pos_phi1.slice(k)).row(0) = pos2_init_mat.col(k).t();
      (IE_neg_phi2.slice(k)).row(0) = neg2_init_mat.col(k).t();
      
      arma::vec pos_part1 = square(IE_list.slice(u) * IE_pos_theta1.slice(k).row(0).t());
      arma::vec pos_part2 = square(IE_list.slice(v) * IE_pos_phi1.slice(k).row(0).t());
      arma::vec neg_part1 = square(IE_list.slice(u) * IE_neg_theta2.slice(k).row(0).t());
      arma::vec neg_part2 = square(IE_list.slice(v) * IE_neg_phi2.slice(k).row(0).t());
      
      arma::vec pos_part = pos_part1 % pos_part2;
      arma::vec neg_part = neg_part1 % neg_part2;
      
      all_interactions.col(k) = pos_part - neg_part;
      
      //Initialize penalty parameter
      
      IE_pen(0,k) = 1;
      
    }
    
    arma::mat accept_MALA(MC, K, fill::zeros);
    accept_MALA.row(0) = vec(K, fill::ones).t();
    
    //Begin MCMC sampling. Remove p_cov part if p_cov = 0.
    
    arma::mat whole_prec(1+p_cov+(p*ME_nspl), 1+p_cov+(p*ME_nspl), fill::zeros);
    whole_prec(0,0) = 0.01;
    whole_prec.submat(1, 1, p_cov, p_cov) = 0.01 * eye(p_cov, p_cov);
    
    for(int m=1; m<MC; ++m){
      
      ////Main Effects and Related Parameters (intercept, sigma^2)
      
      arma::vec R_ME = y - sum(all_interactions, 1);
      
      arma::mat mat_ME_scales(p, p, fill::zeros);
      mat_ME_scales.diag() = ME_scale_stor.row(m-1).t();
      
      arma::mat Psi_ME_inv = kron(mat_ME_scales, SigmaME_inv);
      
      whole_prec.submat(p_cov+1, p_cov+1, 
                        ((p*ME_nspl)+p_cov), ((p*ME_nspl)+p_cov)) = Psi_ME_inv;
      
      arma::vec sampler_ME = maineffects_sampler(R_ME, 
                                                 ME_mat, 
                                                 whole_prec,
                                                 sigmasq_stor(m-1));
      
      alpha_stor(m) = sampler_ME(0);
      cov_effect_stor.row(m) = sampler_ME(span(1, p_cov)).t();
      ME_coeff_stor.row(m) = sampler_ME(span(p_cov+1,p_cov+(p*ME_nspl))).t();
      
      //Interaction Effects and Related Parameters. #####
      
      arma::vec QF13_stor(K, fill::zeros);
      arma::vec QF24_stor(K, fill::zeros);
      arma::vec rej_len_stor(K, fill::zeros);
      
      for(int k=0; k<K; ++k){
        
        if(zero_ind(k) == 1){
          
          // Obtain inverse index (u,v)
          
          int u = map_k_to_uv(k,1);
          int v = map_k_to_uv(k,2);
          
          // Remove index k and sum up other interactions
          
          arma::mat dummy_all_int = all_interactions;
          dummy_all_int.shed_col(k);
          arma::vec all_interaction_sum_notk = sum(dummy_all_int, 1);
          
          ////////// Carry out rejection sampling for interaction k
          
          Rcpp::List rej_pen_obj = rejection_sampler(IE_list.slice(u),
                                                     IE_list.slice(v),
                                                     IE_pen(m-1,k),
                                                     IE_scale_deltasq(m-1),
                                                     IE_scale_tausq1(m-1,k),
                                                     IE_scale_tausq2(m-1,k),
                                                     SigmaInt,
                                                     SigmaInt,
                                                     SigmaInt_inv,
                                                     SigmaInt_inv);
          
          QF13_stor(k) = rej_pen_obj["QF13"];
          QF24_stor(k) = rej_pen_obj["QF24"];
          rej_len_stor(k) = rej_pen_obj["rej_len"];
          
          /////////// Sample parameters for interaction k
          
          arma::vec R_IE_k = y - ((ME_mat * sampler_ME) +
            (all_interaction_sum_notk));
          
          // Define old_param for sampling
          
          arma::vec old_param(n_HMC_param, fill::zeros);
          
          old_param(span(0, IE_nspl - 1)) = IE_pos_theta1.slice(k).row(m-1).t();
          old_param(span(IE_nspl, (2*IE_nspl) - 1)) = IE_neg_theta2.slice(k).row(m-1).t();
          old_param(span(2*IE_nspl, (3*IE_nspl)-1)) = IE_pos_phi1.slice(k).row(m-1).t();
          old_param(span(3*IE_nspl, (4*IE_nspl)-1)) = IE_neg_phi2.slice(k).row(m-1).t();
          
          old_param(4*IE_nspl) = log(IE_pen(m-1,k));
          old_param(4*IE_nspl + 1) = log(IE_scale_tausq1(m-1,k)) / 2.0;
          old_param(4*IE_nspl + 2) = log(IE_scale_tausq2(m-1,k)) / 2.0;
          
          // Do the sampling
          
          arma::mat precond_k_init = pow(c_HMC, 2.0) * eye(n_HMC_param, n_HMC_param);
          arma::mat precond_k_inv_init = eye(n_HMC_param, n_HMC_param) / pow(c_HMC, 2.0);
          
          int len_psi = IE_nspl;
          arma::vec whole_coeff(n_HMC_param, fill::zeros);
          arma::vec whole_coeff_psi(4*len_psi, fill::zeros);
          
          if(m < 5000){
            
            Rcpp::List sq_MALA_k = sq_sampler(R_IE_k, 
                                              IE_list.slice(u), 
                                              IE_list.slice(v), 
                                              SigmaInt_inv,
                                              SigmaInt_inv,
                                              sigmasq_stor(m-1),
                                              IE_scale_deltasq(m-1),
                                              old_param,
                                              eps_MALA(k),
                                              precond_k_init,
                                              precond_k_inv_init,
                                              L_HMC,
                                              rej_pen_obj);
            
            accept_MALA(m,k) = sq_MALA_k["accept"];
            arma::vec sampled_param_k = sq_MALA_k["sampled_param"];
            whole_coeff = sampled_param_k;
            whole_coeff_psi = sampled_param_k(span(0, 4*len_psi - 1));
            
            // Store the sampled gradient 
            
            arma::vec sampled_grad_k = grad_MALA(R_IE_k,
                                                 IE_list.slice(u),
                                                 IE_list.slice(v),
                                                 whole_coeff,
                                                 SigmaInt_inv,
                                                 SigmaInt_inv,
                                                 sigmasq_stor(m-1),
                                                 IE_scale_deltasq(m-1),
                                                 rej_pen_obj);
            
            grad_sample_stor.slice(k).row(m) = sampled_grad_k.t();
            
          }else{
            
            if(m % 500 == 0){
              
              if(precond == 1){
                
                // arma::mat precond_woburnin = grad_sample_stor.slice(k).rows(1000, m-1);
                // 
                // arma::vec precond_mat_mean = mean(precond_woburnin, 0).t();
                // 
                // precond_mat_stor.slice(k) = ((precond_woburnin.t() * precond_woburnin) / m) -
                //   (precond_mat_mean * precond_mat_mean.t());
                
                precond_mat_stor.slice(k) = cov(grad_sample_stor.slice(k).rows(1000, m-1));
                
                precond_mat_inv_stor.slice(k) = inv_sympd(precond_mat_stor.slice(k) +
                  (0.01 * eye(n_HMC_param, n_HMC_param)));
                
              }else{
                
                precond_mat_stor.slice(k) = precond_k_init;
                precond_mat_inv_stor.slice(k) = precond_k_inv_init;
                
              }
              
            }
            
            Rcpp::List sq_MALA_k = sq_sampler(R_IE_k, 
                                              IE_list.slice(u), 
                                              IE_list.slice(v), 
                                              SigmaInt_inv,
                                              SigmaInt_inv,
                                              sigmasq_stor(m-1),
                                              IE_scale_deltasq(m-1),
                                              old_param,
                                              eps_MALA(k),
                                              precond_mat_stor.slice(k),
                                              precond_mat_inv_stor.slice(k),
                                              L_HMC,
                                              rej_pen_obj);
            
            accept_MALA(m,k) = sq_MALA_k["accept"];
            arma::vec sampled_param_k = sq_MALA_k["sampled_param"];
            whole_coeff = sampled_param_k;
            whole_coeff_psi = sampled_param_k(span(0, 4*len_psi - 1));
            
            // Store the sampled gradient 
            
            arma::vec sampled_grad_k = grad_MALA(R_IE_k,
                                                 IE_list.slice(u),
                                                 IE_list.slice(v),
                                                 whole_coeff,
                                                 SigmaInt_inv,
                                                 SigmaInt_inv,
                                                 sigmasq_stor(m-1),
                                                 IE_scale_deltasq(m-1),
                                                 rej_pen_obj);
            
            grad_sample_stor.slice(k).row(m) = sampled_grad_k.t();
            
          }
          
          IE_pos_theta1.slice(k).row(m) = whole_coeff_psi(span(0, len_psi-1)).t();
          IE_neg_theta2.slice(k).row(m) = whole_coeff_psi(span(len_psi, 2*len_psi-1)).t();
          IE_pos_phi1.slice(k).row(m) = whole_coeff_psi(span(2*len_psi, 3*len_psi-1)).t();
          IE_neg_phi2.slice(k).row(m) = whole_coeff_psi(span(3*len_psi, 4*len_psi-1)).t();
          
          IE_pen(m,k) = exp(whole_coeff(4*len_psi));
          IE_scale_tausq1(m,k) = exp(2.0 * whole_coeff(4*len_psi + 1));
          IE_scale_tausq2(m,k) = exp(2.0 * whole_coeff(4*len_psi + 2));
          
          // Update the interaction matrix
          
          arma::vec Ppart1 = square(IE_list.slice(u) * IE_pos_theta1.slice(k).row(m).t());
          arma::vec Ppart2 = square(IE_list.slice(v) * IE_pos_phi1.slice(k).row(m).t());
          arma::vec Npart1 = square(IE_list.slice(u) * IE_neg_theta2.slice(k).row(m).t());
          arma::vec Npart2 = square(IE_list.slice(v) * IE_neg_phi2.slice(k).row(m).t());
          
          arma::vec Ppart = Ppart1 % Ppart2;
          arma::vec Npart = Npart1 % Npart2;
          
          all_interactions.col(k) = Ppart - Npart;
          IE_pen_stor(m,k) = mean(Ppart) * mean(Npart);
          
        }
        
      }
      
      //Sample \lambda_j = scale of \beta_j in ME ####
      
      arma::vec qf_stor(p, fill::ones);
      
      for(int j=0; j<p; ++j){
        
        arma::vec beta_j = ME_coeff_stor(m, span(j*ME_nspl, ((j+1)*ME_nspl)-1)).t();
        arma::mat qf_j = beta_j.t() * (SigmaME_inv * beta_j);
        qf_stor(j) = qf_j(0,0);
        
        ME_scale_stor(m,j) = random_gamma(a_lamb + 
          (0.5*ME_nspl))/(b_lamb + (0.5*qf_stor(j)));
        
        // ME_scale_aux(m,j) = random_gamma(1.0) / (1 + ME_scale_stor(m,j));
        
      }
      
      //// Sample Interaction Effect variance parameters ////
      
      arma::mat qf_stor_int(K, 2, fill::ones);
      
      for(int k=0; k<K; ++k){
        
        if(zero_ind(k) == 1){
          
          arma::mat c_1k1_mat = 0.5*(IE_pos_theta1.slice(k).row(m) * (SigmaInt_inv *
            IE_pos_theta1.slice(k).row(m).t()));
          double c_1k1 = c_1k1_mat(0,0);
          
          arma::mat c_1k2_mat = 0.5*(IE_pos_phi1.slice(k).row(m) * (SigmaInt_inv *
            IE_pos_phi1.slice(k).row(m).t()));
          double c_1k2 = c_1k2_mat(0,0);
          
          double c_1k = c_1k1 + c_1k2;
          
          arma::mat c_2k1_mat = 0.5*(IE_neg_theta2.slice(k).row(m) * (SigmaInt_inv *
            IE_neg_theta2.slice(k).row(m).t()));
          double c_2k1 = c_2k1_mat(0,0);
          
          arma::mat c_2k2_mat = 0.5*(IE_neg_phi2.slice(k).row(m) * (SigmaInt_inv *
            IE_neg_phi2.slice(k).row(m).t()));
          double c_2k2 = c_2k2_mat(0,0);
          
          double c_2k = c_2k1 + c_2k2;
          
          // // Sample \tau_1^2 | a and a | \tau_1^2
          // 
          // double tausq1_rate = (1.0/IE_scale_a(m-1,k)) + 
          //   (c_1k/IE_scale_deltasq(m-1));
          // 
          // IE_scale_tausq1(m,k) = tausq1_rate / random_gamma(IE_nspl + 0.5);
          // IE_scale_a(m,k) = (1 + (1/IE_scale_tausq1(m,k))) / random_gamma(1.0);
          // 
          // // Sample \tau_2^2 | b and b | \tau_2^2
          // 
          // double tausq2_rate = (1.0/IE_scale_b(m-1,k)) + 
          //   (c_2k/IE_scale_deltasq(m-1));
          // 
          // IE_scale_tausq2(m,k) = tausq2_rate / random_gamma(IE_nspl + 0.5);
          // IE_scale_b(m,k) = (1 + (1/IE_scale_tausq2(m,k))) / random_gamma(1.0);
          
          //Store qf/tau^2
          
          qf_stor_int(k,0) = (c_1k + (0.5*QF13_stor(k))) / IE_scale_tausq1(m,k);
          qf_stor_int(k,1) = (c_2k + (0.5*QF24_stor(k))) / IE_scale_tausq2(m,k);
          
        }
        
      }  
      
      // Sample \delta^2 | \nu and \nu | \delta^2
      
      double deltasq_shape = (2*IE_nspl*K) + (2*IE_nspl*accu(rej_len_stor)) + 0.5;
      double deltasq_rate = (accu(qf_stor_int)) + (1/IE_scale_nu(m-1));
      
      IE_scale_deltasq(m) = deltasq_rate / random_gamma(deltasq_shape);
      IE_scale_nu(m) = (1 + (1/IE_scale_deltasq(m))) / random_gamma(1.0);
      
      // Sample error variance
      
      arma::vec R_sigsq = y - ((ME_mat * sampler_ME) + sum(all_interactions,1));
      sigmasq_stor(m) = sigmasq_sampler(R_sigsq, n);
      
      if((m >= 499) && (m < cutoff)){
        
        if(m % 100 == 0){
          
          for(int k1=0; k1<K; ++k1){
            
            if(mean(accept_MALA.rows(m-499,m).col(k1)) <= accept_low){
              
              eps_MALA(k1) = eps_MALA(k1)*accept_scale;
              
            }else if(mean(accept_MALA.rows(m-499,m).col(k1)) >= accept_high){
              
              eps_MALA(k1) = eps_MALA(k1)/accept_scale;
              
            }else{
              
              eps_MALA(k1) = eps_MALA(k1);
              
            }
            
          }
          
        }
        
      }  
      
      if((m+1) % 1000 == 0){
        
        Rcpp::Rcout << "Monte Carlo Sample: " << m+1 << std::endl;
        
      }
      
    }  
    
    return Rcpp::List::create(Rcpp::Named("intercept") = alpha_stor,
                              Rcpp::Named("error_var") = sigmasq_stor,
                              Rcpp::Named("ME_coeff") = ME_coeff_stor,
                              Rcpp::Named("ME_scale") = ME_scale_stor,
                              Rcpp::Named("IE_pos_coeff_part1") = IE_pos_theta1,
                              Rcpp::Named("IE_pos_coeff_part2") = IE_pos_phi1,
                              Rcpp::Named("IE_neg_coeff_part1") = IE_neg_theta2,
                              Rcpp::Named("IE_neg_coeff_part2") =  IE_neg_phi2,
                              Rcpp::Named("cov_effects") = cov_effect_stor,
                              Rcpp::Named("ME_mat") = ME_mat,
                              Rcpp::Named("IE_list") = IE_list,
                              Rcpp::Named("Accept_Prop") = accept_MALA,
                              Rcpp::Named("HMC_epsilon") = eps_MALA,
                              Rcpp::Named("IE_penalty") = IE_pen,
                              Rcpp::Named("IE_penalty_stor") = IE_pen_stor,
                              Rcpp::Named("IE_deltasq") = IE_scale_deltasq,
                              Rcpp::Named("IE_tausq1") = IE_scale_tausq1,
                              Rcpp::Named("IE_tausq2") = IE_scale_tausq2,
                              Rcpp::Named("precond_mat_stor") = precond_mat_stor);
  }else{
    
    //// Initialize ////
    
    // Intercept and sigma^2
    
    alpha_stor(0) = init_values["intercept_est"];
    sigmasq_stor(0) = init_values["sigmasq_est"];
    
    // Main effects
    
    arma::vec ME_coeff_init_vec = init_values["ME_coeff_est"];
    
    ME_coeff_stor.row(0) = ME_coeff_init_vec.t();
    
    // Interaction effects
    
    arma::mat pos1_init_mat = init_values["IE_pos_part1_coeff_est"];
    arma::mat pos2_init_mat = init_values["IE_pos_part2_coeff_est"];
    arma::mat neg1_init_mat = init_values["IE_neg_part1_coeff_est"];
    arma::mat neg2_init_mat = init_values["IE_neg_part2_coeff_est"];
    
    IE_scale_deltasq(0) = 0.01;
    
    for(int k=0; k<K; ++k){
      
      // Obtain inverse index (u,v)
      
      int u = map_k_to_uv(k,1);
      int v = map_k_to_uv(k,2);
      
      (IE_pos_theta1.slice(k)).row(0) = pos1_init_mat.col(k).t();
      (IE_neg_theta2.slice(k)).row(0) = neg1_init_mat.col(k).t();
      (IE_pos_phi1.slice(k)).row(0) = pos2_init_mat.col(k).t();
      (IE_neg_phi2.slice(k)).row(0) = neg2_init_mat.col(k).t();
      
      arma::vec pos_part1 = square(IE_list.slice(u) * IE_pos_theta1.slice(k).row(0).t());
      arma::vec pos_part2 = square(IE_list.slice(v) * IE_pos_phi1.slice(k).row(0).t());
      arma::vec neg_part1 = square(IE_list.slice(u) * IE_neg_theta2.slice(k).row(0).t());
      arma::vec neg_part2 = square(IE_list.slice(v) * IE_neg_phi2.slice(k).row(0).t());
      
      arma::vec pos_part = pos_part1 % pos_part2;
      arma::vec neg_part = neg_part1 % neg_part2;
      
      all_interactions.col(k) = pos_part - neg_part;
      
      //Initialize penalty parameter
      
      IE_pen(0,k) = 1;
      
    }
    
    arma::mat accept_MALA(MC, K, fill::zeros);
    accept_MALA.row(0) = vec(K, fill::ones).t();
    
    //Begin MCMC sampling. Remove p_cov part if p_cov = 0.
    
    arma::mat whole_prec(1+p_cov+(p*ME_nspl), 1+p_cov+(p*ME_nspl), fill::zeros);
    whole_prec(0,0) = 0.01;
    
    for(int m=1; m<MC; ++m){
      
      ////Main Effects and Related Parameters (intercept, sigma^2)
      
      arma::vec R_ME = y - sum(all_interactions, 1);
      
      arma::mat mat_ME_scales(p, p, fill::zeros);
      mat_ME_scales.diag() = ME_scale_stor.row(m-1).t();
      
      arma::mat Psi_ME_inv = kron(mat_ME_scales, SigmaME_inv);
      
      whole_prec.submat(p_cov+1, p_cov+1, 
                        ((p*ME_nspl)+p_cov), ((p*ME_nspl)+p_cov)) = Psi_ME_inv;
      
      arma::vec sampler_ME = maineffects_sampler(R_ME, 
                                                 ME_mat, 
                                                 whole_prec,
                                                 sigmasq_stor(m-1));
      
      alpha_stor(m) = sampler_ME(0);
      ME_coeff_stor.row(m) = sampler_ME(span(p_cov+1,p_cov+(p*ME_nspl))).t();
      
      //Interaction Effects and Related Parameters. #####
      
      arma::vec QF13_stor(K, fill::zeros);
      arma::vec QF24_stor(K, fill::zeros);
      arma::vec rej_len_stor(K, fill::zeros);
      
      for(int k=0; k<K; ++k){
        
        if(zero_ind(k) == 1){
          
          // Obtain inverse index (u,v)
          
          int u = map_k_to_uv(k,1);
          int v = map_k_to_uv(k,2);
          
          // Remove index k and sum up other interactions
          
          arma::mat dummy_all_int = all_interactions;
          dummy_all_int.shed_col(k);
          arma::vec all_interaction_sum_notk = sum(dummy_all_int, 1);
          
          ////////// Carry out rejection sampling for interaction k
          
          Rcpp::List rej_pen_obj = rejection_sampler(IE_list.slice(u),
                                                     IE_list.slice(v),
                                                     IE_pen(m-1,k),
                                                     IE_scale_deltasq(m-1),
                                                     IE_scale_tausq1(m-1,k),
                                                     IE_scale_tausq2(m-1,k),
                                                     SigmaInt,
                                                     SigmaInt,
                                                     SigmaInt_inv,
                                                     SigmaInt_inv);
          
          QF13_stor(k) = rej_pen_obj["QF13"];
          QF24_stor(k) = rej_pen_obj["QF24"];
          rej_len_stor(k) = rej_pen_obj["rej_len"];
          
          /////////// Sample parameters for interaction k
          
          arma::vec R_IE_k = y - ((ME_mat * sampler_ME) + (all_interaction_sum_notk));
          
          // Define old_param for sampling
          
          arma::vec old_param(n_HMC_param, fill::zeros);
          
          old_param(span(0, IE_nspl - 1)) = IE_pos_theta1.slice(k).row(m-1).t();
          old_param(span(IE_nspl, (2*IE_nspl) - 1)) = IE_neg_theta2.slice(k).row(m-1).t();
          old_param(span(2*IE_nspl, (3*IE_nspl)-1)) = IE_pos_phi1.slice(k).row(m-1).t();
          old_param(span(3*IE_nspl, (4*IE_nspl)-1)) = IE_neg_phi2.slice(k).row(m-1).t();
          
          old_param(4*IE_nspl) = log(IE_pen(m-1,k));
          old_param(4*IE_nspl + 1) = log(IE_scale_tausq1(m-1,k)) / 2.0;
          old_param(4*IE_nspl + 2) = log(IE_scale_tausq2(m-1,k)) / 2.0;
          
          // Do the sampling
          
          arma::mat precond_k_init = pow(c_HMC, 2.0) * eye(n_HMC_param, n_HMC_param);
          arma::mat precond_k_inv_init = eye(n_HMC_param, n_HMC_param) / pow(c_HMC, 2.0);
          
          int len_psi = IE_nspl;
          arma::vec whole_coeff(n_HMC_param, fill::zeros);
          arma::vec whole_coeff_psi(4*len_psi, fill::zeros);
          
          if(m < 5000){
            
            Rcpp::List sq_MALA_k = sq_sampler(R_IE_k, 
                                              IE_list.slice(u), 
                                              IE_list.slice(v), 
                                              SigmaInt_inv,
                                              SigmaInt_inv,
                                              sigmasq_stor(m-1),
                                              IE_scale_deltasq(m-1),
                                              old_param,
                                              eps_MALA(k),
                                              precond_k_init,
                                              precond_k_inv_init,
                                              L_HMC,
                                              rej_pen_obj);
            
            accept_MALA(m,k) = sq_MALA_k["accept"];
            arma::vec sampled_param_k = sq_MALA_k["sampled_param"];
            whole_coeff = sampled_param_k;
            whole_coeff_psi = sampled_param_k(span(0, 4*len_psi - 1));
            
            // Store the sampled gradient 
            
            arma::vec sampled_grad_k = grad_MALA(R_IE_k,
                                                 IE_list.slice(u),
                                                 IE_list.slice(v),
                                                 whole_coeff,
                                                 SigmaInt_inv,
                                                 SigmaInt_inv,
                                                 sigmasq_stor(m-1),
                                                 IE_scale_deltasq(m-1),
                                                 rej_pen_obj);
            
            grad_sample_stor.slice(k).row(m) = sampled_grad_k.t();
            
          }else{
            
            if(m % 500 == 0){
              
              if(precond == 1){
                
                // arma::mat precond_woburnin = grad_sample_stor.slice(k).rows(1000, m-1);
                // 
                // arma::vec precond_mat_mean = mean(precond_woburnin, 0).t();
                // 
                // precond_mat_stor.slice(k) = ((precond_woburnin.t() * precond_woburnin) / m) -
                //   (precond_mat_mean * precond_mat_mean.t());
                
                precond_mat_stor.slice(k) = cov(grad_sample_stor.slice(k).rows(1000, m-1));
                
                precond_mat_inv_stor.slice(k) = inv_sympd(precond_mat_stor.slice(k) +
                  (0.01 * eye(n_HMC_param, n_HMC_param)));
                
              }else{
                
                precond_mat_stor.slice(k) = precond_k_init;
                precond_mat_inv_stor.slice(k) = precond_k_inv_init;
                
              }
              
            }
            
            Rcpp::List sq_MALA_k = sq_sampler(R_IE_k, 
                                              IE_list.slice(u), 
                                              IE_list.slice(v), 
                                              SigmaInt_inv,
                                              SigmaInt_inv,
                                              sigmasq_stor(m-1),
                                              IE_scale_deltasq(m-1),
                                              old_param,
                                              eps_MALA(k),
                                              precond_mat_stor.slice(k),
                                              precond_mat_inv_stor.slice(k),
                                              L_HMC,
                                              rej_pen_obj);
            
            accept_MALA(m,k) = sq_MALA_k["accept"];
            arma::vec sampled_param_k = sq_MALA_k["sampled_param"];
            whole_coeff = sampled_param_k;
            whole_coeff_psi = sampled_param_k(span(0, 4*len_psi - 1));
            
            // Store the sampled gradient 
            
            arma::vec sampled_grad_k = grad_MALA(R_IE_k,
                                                 IE_list.slice(u),
                                                 IE_list.slice(v),
                                                 whole_coeff,
                                                 SigmaInt_inv,
                                                 SigmaInt_inv,
                                                 sigmasq_stor(m-1),
                                                 IE_scale_deltasq(m-1),
                                                 rej_pen_obj);
            
            grad_sample_stor.slice(k).row(m) = sampled_grad_k.t();
            
          }
          
          IE_pos_theta1.slice(k).row(m) = whole_coeff_psi(span(0, len_psi-1)).t();
          IE_neg_theta2.slice(k).row(m) = whole_coeff_psi(span(len_psi, 2*len_psi-1)).t();
          IE_pos_phi1.slice(k).row(m) = whole_coeff_psi(span(2*len_psi, 3*len_psi-1)).t();
          IE_neg_phi2.slice(k).row(m) = whole_coeff_psi(span(3*len_psi, 4*len_psi-1)).t();
          
          IE_pen(m,k) = exp(whole_coeff(4*len_psi));
          IE_scale_tausq1(m,k) = exp(2.0 * whole_coeff(4*len_psi + 1));
          IE_scale_tausq2(m,k) = exp(2.0 * whole_coeff(4*len_psi + 2));
          
          // Update the interaction matrix
          
          arma::vec Ppart1 = square(IE_list.slice(u) * IE_pos_theta1.slice(k).row(m).t());
          arma::vec Ppart2 = square(IE_list.slice(v) * IE_pos_phi1.slice(k).row(m).t());
          arma::vec Npart1 = square(IE_list.slice(u) * IE_neg_theta2.slice(k).row(m).t());
          arma::vec Npart2 = square(IE_list.slice(v) * IE_neg_phi2.slice(k).row(m).t());
          
          arma::vec Ppart = Ppart1 % Ppart2;
          arma::vec Npart = Npart1 % Npart2;
          
          all_interactions.col(k) = Ppart - Npart; 
          IE_pen_stor(m,k) = mean(Ppart) * mean(Npart);
          
        }
        
      }
      
      // Rcpp::Rcout << rej_len_stor << std::endl;
      
      //Sample \lambda_j = scale of \beta_j in ME ####
      
      arma::vec qf_stor(p, fill::ones);
      
      for(int j=0; j<p; ++j){
        
        arma::vec beta_j = ME_coeff_stor(m, span(j*ME_nspl, ((j+1)*ME_nspl)-1)).t();
        arma::mat qf_j = beta_j.t() * (SigmaME_inv * beta_j);
        qf_stor(j) = qf_j(0,0);
        
        ME_scale_stor(m,j) = random_gamma(a_lamb + 
          (0.5*ME_nspl))/(b_lamb + (0.5*qf_stor(j)));
        
        // ME_scale_aux(m,j) = random_gamma(1.0) / (1 + ME_scale_stor(m,j));
        
      }
      
      //// Sample Interaction Effect variance parameters ////
      
      arma::mat qf_stor_int(K, 2, fill::ones);
      
      for(int k=0; k<K; ++k){
        
        if(zero_ind(k) == 1){
          
          arma::mat c_1k1_mat = 0.5*(IE_pos_theta1.slice(k).row(m) * (SigmaInt_inv *
            IE_pos_theta1.slice(k).row(m).t()));
          double c_1k1 = c_1k1_mat(0,0);
          
          arma::mat c_1k2_mat = 0.5*(IE_pos_phi1.slice(k).row(m) * (SigmaInt_inv *
            IE_pos_phi1.slice(k).row(m).t()));
          double c_1k2 = c_1k2_mat(0,0);
          
          double c_1k = c_1k1 + c_1k2;
          
          arma::mat c_2k1_mat = 0.5*(IE_neg_theta2.slice(k).row(m) * (SigmaInt_inv *
            IE_neg_theta2.slice(k).row(m).t()));
          double c_2k1 = c_2k1_mat(0,0);
          
          arma::mat c_2k2_mat = 0.5*(IE_neg_phi2.slice(k).row(m) * (SigmaInt_inv *
            IE_neg_phi2.slice(k).row(m).t()));
          double c_2k2 = c_2k2_mat(0,0);
          
          double c_2k = c_2k1 + c_2k2;
          
          // // Sample \tau_1^2 | a and a | \tau_1^2
          // 
          // double tausq1_rate = (1.0/IE_scale_a(m-1,k)) + 
          //   (c_1k/IE_scale_deltasq(m-1));
          // 
          // IE_scale_tausq1(m,k) = tausq1_rate / random_gamma(IE_nspl + 0.5);
          // IE_scale_a(m,k) = (1 + (1/IE_scale_tausq1(m,k))) / random_gamma(1.0);
          // 
          // // Sample \tau_2^2 | b and b | \tau_2^2
          // 
          // double tausq2_rate = (1.0/IE_scale_b(m-1,k)) + 
          //   (c_2k/IE_scale_deltasq(m-1));
          // 
          // IE_scale_tausq2(m,k) = tausq2_rate / random_gamma(IE_nspl + 0.5);
          // IE_scale_b(m,k) = (1 + (1/IE_scale_tausq2(m,k))) / random_gamma(1.0);
          
          //Store qf/tau^2
          
          qf_stor_int(k,0) = (c_1k + (0.5*QF13_stor(k))) / IE_scale_tausq1(m,k);
          qf_stor_int(k,1) = (c_2k + (0.5*QF24_stor(k))) / IE_scale_tausq2(m,k);
          
        }
        
      }  
      
      // Sample \delta^2 | \nu and \nu | \delta^2
      
      double deltasq_shape = (2*IE_nspl*K) + (2*IE_nspl*accu(rej_len_stor)) + 0.5;
      double deltasq_rate = (accu(qf_stor_int)) + (1/IE_scale_nu(m-1));
      
      IE_scale_deltasq(m) = deltasq_rate / random_gamma(deltasq_shape);
      //IE_scale_deltasq(m) = 1.0;
      IE_scale_nu(m) = (1 + (1/IE_scale_deltasq(m))) / random_gamma(1.0);
      
      // Sample error variance
      
      arma::vec R_sigsq = y - ((ME_mat * sampler_ME) + sum(all_interactions,1));
      sigmasq_stor(m) = sigmasq_sampler(R_sigsq, n);
      
      if((m >= 499) && (m < cutoff)){
        
        if(m % 100 == 0){
          
          for(int k1=0; k1<K; ++k1){
            
            if(mean(accept_MALA.rows(m-499,m).col(k1)) <= accept_low){
              
              eps_MALA(k1) = eps_MALA(k1)*accept_scale;
              
            }else if(mean(accept_MALA.rows(m-499,m).col(k1)) >= accept_high){
              
              eps_MALA(k1) = eps_MALA(k1)/accept_scale;
              
            }else{
              
              eps_MALA(k1) = eps_MALA(k1);
              
            }
            
          }
          
        }
        
      }  
      
      if((m+1) % 1000 == 0){
        
        Rcpp::Rcout << "Monte Carlo Sample: " << m+1 << std::endl;
        
      }
      
    }  
    
    return Rcpp::List::create(Rcpp::Named("intercept") = alpha_stor,
                              Rcpp::Named("error_var") = sigmasq_stor,
                              Rcpp::Named("ME_coeff") = ME_coeff_stor,
                              Rcpp::Named("ME_scale") = ME_scale_stor,
                              Rcpp::Named("IE_pos_coeff_part1") = IE_pos_theta1,
                              Rcpp::Named("IE_pos_coeff_part2") = IE_pos_phi1,
                              Rcpp::Named("IE_neg_coeff_part1") = IE_neg_theta2,
                              Rcpp::Named("IE_neg_coeff_part2") =  IE_neg_phi2,
                              Rcpp::Named("ME_mat") = ME_mat,
                              Rcpp::Named("IE_list") = IE_list,
                              Rcpp::Named("Accept_Prop") = accept_MALA,
                              Rcpp::Named("HMC_epsilon") = eps_MALA,
                              Rcpp::Named("IE_penalty") = IE_pen,
                              Rcpp::Named("IE_penalty_stor") = IE_pen_stor,
                              Rcpp::Named("IE_deltasq") = IE_scale_deltasq,
                              Rcpp::Named("IE_tausq1") = IE_scale_tausq1,
                              Rcpp::Named("IE_tausq2") = IE_scale_tausq2,
                              Rcpp::Named("precond_mat_stor") = precond_mat_stor);
    
  }
  
}
