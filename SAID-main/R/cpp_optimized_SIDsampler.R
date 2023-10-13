### Univariate Penalty

penmatt<-function(M)
{

  Amat = matrix(0, nrow = M, ncol = M)

  Amat[1,] = c(1, rep(0, M-1))

  for(i in 2:M)
  {

    Amat[i,] = c(rep(0, i-2), -1, 1, rep(0, M-i))

  }

  return(t(Amat) %*% Amat)

}

### Bivariate Penalty

# penmatt<-function(M)
# {
# 
#   Amat = matrix(0, nrow = M, ncol = M)
# 
#   Amat[1,] = c(1, rep(0, M-1))
#   Amat[2,] = c(-2, 1, rep(0, M-2))
# 
#   for(i in 3:M)
#   {
# 
#     Amat[i,] = c(rep(0, i-3), 1, -2, 1, rep(0, M-i))
# 
#   }
# 
#   return(t(Amat) %*% Amat)
# 
# }

SIMinit<-function(y, X, Z = NULL, ME_list, IE_list, zero_ind, nugget,
                  me_integral_constraint = TRUE)
{
  
  library(MASS)
  library(randomForest)
  
  #### Fit a linear model ####
  
  p_cov = ifelse(is.null(Z) == TRUE, 0, dim(Z)[2])
  
  if(p_cov > 0)
  {
    
    lm_mod = lm(y ~ Z - 1)
    cov_effect_est = lm_mod$coefficients
    
    R = y - (Z %*% cov_effect_est)
    
  }else
  {
    
    R = y
    
  }
  
  #### Fit random forest on R ####
  
  RFmodel = randomForest::randomForest(x = X, y = R)
  
  if(me_integral_constraint == TRUE)
  {
    
    ## Extract the main effect coefficients ##
    
    n = dim(X)[1]
    p = dim(X)[2]
    nspl_ME = dim(ME_list)[2]
    
    overall_int = rep(0, p)
    
    ME_coeff_est_mat = matrix(0, nrow = nspl_ME, ncol = p)
    
    for(j in 1:p)
    {
      
      Xtest_j = matrix(0, nrow = n, ncol = p)
      Xtest_j[,j] = X[,j]
      RF_ME_j_11 = predict(RFmodel, newdata = Xtest_j)
      
      overall_int[j] = mean(RF_ME_j_11)
      
      ME_j = RF_ME_j_11 - overall_int[j]
      
      M_j = ME_list[,,j]
      
      # ME_coeff_est_mat[,j] = ginv(M_j) %*% ME_j
        
      ME_coeff_est_mat[,j] = solve((t(M_j) %*% M_j) + (nugget * diag(nspl_ME))) %*%
        t(M_j) %*% ME_j
      
    }
    
    ## Extract the intercept term ##
    
    RF_intercept = predict(RFmodel, newdata = rep(0, p))
    intercept_est = sum(overall_int) - ((p-1) * RF_intercept)
    
    ## Extract the interaction effect coefficients ##
    
    nspl_IE = dim(IE_list)[2]
    pos_coeff_part1_int = matrix(rnorm(nspl_IE * choose(p,2)), nrow = nspl_IE, ncol = choose(p,2))
    pos_coeff_part2_int = matrix(rnorm(nspl_IE * choose(p,2)), nrow = nspl_IE, ncol = choose(p,2))
    neg_coeff_part1_int = matrix(rnorm(nspl_IE * choose(p,2)), nrow = nspl_IE, ncol = choose(p,2))
    neg_coeff_part2_int = matrix(rnorm(nspl_IE * choose(p,2)), nrow = nspl_IE, ncol = choose(p,2))
    
    for(k in 1:choose(p,2))
    {
     
      if(zero_ind[k] == 0)
      {
       
        pos_coeff_part1_int[,k] = rep(0, nspl_IE)
        pos_coeff_part2_int[,k] = rep(0, nspl_IE)
        neg_coeff_part1_int[,k] = rep(0, nspl_IE)
        neg_coeff_part2_int[,k] = rep(0, nspl_IE)
        
      }
      
    }  
    
    # for(k in 1:choose(p,2))
    # {
    # 
    #   if(zero_ind[k] == 1)
    #   {
    # 
    #     #Obtain inverse index (u,v)
    # 
    #     quad_dis = (2*p - 1)^2 - 8*k
    #     u = ceiling(0.5*((2*p-1) - quad_dis^0.5))
    #     v = p + k - (u*(p - 0.5*(u+1)))
    # 
    #     # Extract (u,v)-th interaction estimate #
    # 
    #     S_u = IE_list[,,u]
    #     S_v = IE_list[,,v]
    # 
    #     mat_uv_11 = matrix(0, nrow = n, ncol = p)
    #     mat_uv_10 = matrix(0, nrow = n, ncol = p)
    #     mat_uv_01 = matrix(0, nrow = n, ncol = p)
    # 
    #     mat_uv_11[,u] = X[,u]
    #     mat_uv_11[,v] = X[,v]
    # 
    #     mat_uv_10[,u] = X[,u]
    #     mat_uv_01[,v] = X[,v]
    # 
    #     int_est_uv_11 = predict(RFmodel, newdata = mat_uv_11)
    #     int_est_uv_10 = predict(RFmodel, newdata = mat_uv_10)
    #     int_est_uv_01 = predict(RFmodel, newdata = mat_uv_01)
    # 
    #     int_est_uv = int_est_uv_11 -
    #       int_est_uv_10 -
    #       int_est_uv_01 +
    #       rep(RF_intercept, n)
    # 
    #     # Extract positive part coefficients #
    # 
    #     pos_uv = sapply(int_est_uv, max, 0)
    # 
    #     pos_uv_part1 = pos_uv^(0.25)
    #     pos_uv_part2 = pos_uv^(0.25)
    # 
    #     # pos_coeff_part1_int[,k] = solve((t(S_u) %*% S_u) + (nugget * diag(nspl_IE))) %*%
    #     #   t(S_u) %*% pos_uv_part1
    #     #
    #     # pos_coeff_part2_int[,k] = solve((t(S_v) %*% S_v) + (nugget * diag(nspl_IE))) %*%
    #     #   t(S_v) %*% pos_uv_part2
    # 
    #     pos_coeff_part1_int[,k] = ginv(S_u) %*% pos_uv_part1
    # 
    #     pos_coeff_part2_int[,k] = ginv(S_v) %*% pos_uv_part2
    # 
    #     # Extract negative part coefficients #
    # 
    #     neg_uv = sapply(-int_est_uv, max, 0)
    # 
    #     neg_uv_part1 = neg_uv^(0.25)
    #     neg_uv_part2 = neg_uv^(0.25)
    # 
    #     # neg_coeff_part1_int[,k] = solve((t(S_u) %*% S_u) + (nugget * diag(nspl_IE))) %*%
    #     #   t(S_u) %*% neg_uv_part1
    #     #
    #     # neg_coeff_part2_int[,k] = solve((t(S_v) %*% S_v) + (nugget * diag(nspl_IE))) %*%
    #     #   t(S_v) %*% neg_uv_part2
    # 
    #     neg_coeff_part1_int[,k] = ginv(S_u) %*% neg_uv_part1
    # 
    #     neg_coeff_part2_int[,k] = ginv(S_v) %*% neg_uv_part2
    # 
    #   }
    # 
    # }
    
    ## Extract the variance estimate ##
    
    whole_sur_est = predict(RFmodel, newdata = X)
    
    sigmasq_est = mean((y - whole_sur_est)^2)
    
    if(p_cov > 0)
    {
      
      init_list = list("cov_effect_est" = cov_effect_est,
                       "intercept_est" = intercept_est,
                       "sigmasq_est" = sigmasq_est,
                       "ME_coeff_est" = c(ME_coeff_est_mat),
                       "IE_pos_part1_coeff_est" = pos_coeff_part1_int,
                       "IE_pos_part2_coeff_est" = pos_coeff_part2_int,
                       "IE_neg_part1_coeff_est" = neg_coeff_part1_int,
                       "IE_neg_part2_coeff_est" = neg_coeff_part2_int)
      
    }else
    {
      
      init_list = list("intercept_est" = intercept_est,
                       "sigmasq_est" = sigmasq_est,
                       "ME_coeff_est" = c(ME_coeff_est_mat),
                       "IE_pos_part1_coeff_est" = pos_coeff_part1_int,
                       "IE_pos_part2_coeff_est" = pos_coeff_part2_int,
                       "IE_neg_part1_coeff_est" = neg_coeff_part1_int,
                       "IE_neg_part2_coeff_est" = neg_coeff_part2_int)
      
    }
    
    return(init_list)
    
  }else
  {
    
    ## Extract the main effect coefficients ##
    
    n = dim(X)[1]
    p = dim(X)[2]
    
    ## Extract the intercept ##
    
    RF_intercept = predict(RFmodel, newdata = rep(0, p))
    intercept_est = RF_intercept
    
    nspl_ME = dim(ME_list)[2]
    
    ME_coeff_est_mat = matrix(0, nrow = nspl_ME, ncol = p)
    
    for(j in 1:p)
    {
      
      Xtest_j = matrix(0, nrow = n, ncol = p)
      Xtest_j[,j] = X[,j]
      RF_ME_j_11 = predict(RFmodel, newdata = Xtest_j)
      
      ME_j = RF_ME_j_11 - rep(intercept_est, n)
      
      M_j = ME_list[,,j]
      
      ME_coeff_est_mat[,j] = solve((t(M_j) %*% M_j) + (nugget * diag(nspl_ME))) %*%
        t(M_j) %*% ME_j
      
    }
    
    ## Extract the interaction effect coefficients ##
    
    nspl_IE = dim(IE_list)[2]
    pos_coeff_part1_int = matrix(rnorm(nspl_IE * choose(p,2)), nrow = nspl_IE, ncol = choose(p,2))
    pos_coeff_part2_int = matrix(rnorm(nspl_IE * choose(p,2)), nrow = nspl_IE, ncol = choose(p,2))
    neg_coeff_part1_int = matrix(rnorm(nspl_IE * choose(p,2)), nrow = nspl_IE, ncol = choose(p,2))
    neg_coeff_part2_int = matrix(rnorm(nspl_IE * choose(p,2)), nrow = nspl_IE, ncol = choose(p,2))
    
    for(k in 1:choose(p,2))
    {
     
      if(zero_ind[k] == 0)
      {
       
        pos_coeff_part1_int[,k] = rep(0, nspl_IE)
        pos_coeff_part2_int[,k] = rep(0, nspl_IE)
        neg_coeff_part1_int[,k] = rep(0, nspl_IE)
        neg_coeff_part2_int[,k] = rep(0, nspl_IE)
        
      }
      
    }  
    
#     for(k in 1:choose(p,2))
#     {
      
#       if(zero_ind[k] == 1)
#       {
        
#         #Obtain inverse index (u,v)
        
#         quad_dis = (2*p - 1)^2 - 8*k
#         u = ceiling(0.5*((2*p-1) - quad_dis^0.5))
#         v = p + k - (u*(p - 0.5*(u+1)))
        
#         # Extract (u,v)-th interaction estimate #
        
#         S_u = IE_list[,,u]
#         S_v = IE_list[,,v]
        
#         mat_uv_11 = matrix(0, nrow = n, ncol = p)
#         mat_uv_10 = matrix(0, nrow = n, ncol = p)
#         mat_uv_01 = matrix(0, nrow = n, ncol = p)
        
#         mat_uv_11[,u] = X[,u]
#         mat_uv_11[,v] = X[,v]
        
#         mat_uv_10[,u] = X[,u]
#         mat_uv_01[,v] = X[,v]
        
#         int_est_uv_11 = predict(RFmodel, newdata = mat_uv_11)
#         int_est_uv_10 = predict(RFmodel, newdata = mat_uv_10)
#         int_est_uv_01 = predict(RFmodel, newdata = mat_uv_01)
        
#         int_est_uv = int_est_uv_11 -
#           int_est_uv_10 -
#           int_est_uv_01 +
#           rep(intercept_est, n)
        
#         # Extract positive part coefficients #
        
#         pos_uv = sapply(int_est_uv, max, 0)
        
#         pos_uv_part1 = pos_uv^(0.25)
#         pos_uv_part2 = pos_uv^(0.25)
        
#         pos_coeff_part1_int[,k] = solve((t(S_u) %*% S_u) + (nugget * diag(nspl_IE))) %*%
#           t(S_u) %*% pos_uv_part1
        
#         pos_coeff_part2_int[,k] = solve((t(S_v) %*% S_v) + (nugget * diag(nspl_IE))) %*%
#           t(S_v) %*% pos_uv_part2
        
#         # Extract negative part coefficients #
        
#         neg_uv = sapply(-int_est_uv, max, 0)
        
#         neg_uv_part1 = neg_uv^(0.25)
#         neg_uv_part2 = neg_uv^(0.25)
        
#         neg_coeff_part1_int[,k] = solve((t(S_u) %*% S_u) + (nugget * diag(nspl_IE))) %*%
#           t(S_u) %*% neg_uv_part1
        
#         neg_coeff_part2_int[,k] = solve((t(S_v) %*% S_v) + (nugget * diag(nspl_IE))) %*%
#           t(S_v) %*% neg_uv_part2
        
#       }
      
#     }
    
    ## Extract the variance estimate ##
    
    whole_sur_est = predict(RFmodel, newdata = X)
    
    sigmasq_est = mean((y - whole_sur_est)^2)
    
    if(p_cov > 0)
    {
      
      init_list = list("cov_effect_est" = cov_effect_est,
                       "intercept_est" = intercept_est,
                       "sigmasq_est" = sigmasq_est,
                       "ME_coeff_est" = c(ME_coeff_est_mat),
                       "IE_pos_part1_coeff_est" = pos_coeff_part1_int,
                       "IE_pos_part2_coeff_est" = pos_coeff_part2_int,
                       "IE_neg_part1_coeff_est" = neg_coeff_part1_int,
                       "IE_neg_part2_coeff_est" = neg_coeff_part2_int)
      
    }else
    {
      
      init_list = list("intercept_est" = intercept_est,
                       "sigmasq_est" = sigmasq_est,
                       "ME_coeff_est" = c(ME_coeff_est_mat),
                       "IE_pos_part1_coeff_est" = pos_coeff_part1_int,
                       "IE_pos_part2_coeff_est" = pos_coeff_part2_int,
                       "IE_neg_part1_coeff_est" = neg_coeff_part1_int,
                       "IE_neg_part2_coeff_est" = neg_coeff_part2_int)
      
    }
    
    return(init_list)
    
  }
  
  
  
}

##M+4 = #splines for main effects, N+3 = #splines for interaction tensor products

SAIDsampler<-function(y,
                     X, 
                     Z = NULL, 
                     K_ME = 5,
                     K_IE = 2, 
                     a_lamb = 0.001,
                     b_lamb = 0.001,
                     eps_MALA = rep(0.1, choose(dim(X)[2], 2)),
                     c_HMC = 1, 
                     L_HMC = 5, 
                     MC = 10000,
                     zero_ind = rep(1, choose(dim(X)[2], 2)),
                     me_integral_constraint = TRUE,
                     cutoff = 0.5*MC,
                     accept_low = 0.50,
                     accept_high = 0.60,
                     accept_scale = 0.8,
                     precond = 0){
  
  library(MASS)
  library(splines)
  library(splines2)
  library(mvtnorm)
  
  #### Samples from the model
  #### E(y_i | X_{1:p}) = \alpha +
  #### \sum_{j=1}^{p} f_j(X_{j,i}) +
  #### \sum_{u < v} h_{uv}(X_{u, i}, X_{v, i}) +
  #### \epsilon_i.
  
  n = dim(X)[1]
  p = dim(X)[2]
  
  if(is.null(Z) == TRUE)
  {
    
    p_cov = 0
    
  }else
  {
    
    p_cov = dim(Z)[2]
    
  }
  
  #### B-Splines for computation ####
  
  if(me_integral_constraint == TRUE)
  {
    
    nspl_ME = K_ME + 4
    
  }else
  {
    
    nspl_ME = K_ME + 3
    
  }
  
  nspl_IE = K_IE + 3
  
  ME_list = array(0, dim = c(n, nspl_ME, p))
  IE_list = array(0, dim = c(n, nspl_IE, p))
  ME_subtract = matrix(0, nrow = p, ncol = nspl_ME)
  
  S_ME_inv = array(0, dim = c(nspl_ME, nspl_ME, p))
  S_IE_inv = array(0, dim = c(nspl_IE, nspl_IE, p))
  
  ME_knots_stor = matrix(0, nrow = p, ncol = K_ME)
  IE_knots_stor = matrix(0, nrow = p, ncol = K_IE)
  
  ind = 1
  
  for(ind in 1:p)
  {
    
    ### For main effects ###
    
    #me_knots = seq(1/M, 1-(1/M), length.out = M)
    
    quantile_seq_ME = seq(0, 1, by = 1/(K_ME+1))
    quantile_seq_ME = quantile_seq_ME[-c(1,K_ME+2)]
    
    me_knots = quantile(X[,ind], quantile_seq_ME)
    
    ME_knots_stor[ind, ] = me_knots
    
    if(me_integral_constraint == TRUE)
    {
    
      me_spl = bSpline(x = X[,ind], knots = me_knots, intercept = TRUE)
      ME_subtract[ind,] = colMeans(me_spl)
      final_Xmat_ME = sweep(me_spl, 2,  ME_subtract[ind,])
    
    }else
    {
    
      me_spl = bSpline(x = X[,ind], knots = me_knots, intercept = FALSE)
      final_Xmat_ME = me_spl
      
    }
    
    ME_list[,,ind] = final_Xmat_ME
    
    S_ME_inv[,,ind] = penmatt(nspl_ME)
    
    ### For interaction effects ###
    
    quantile_seq_IE = seq(0, 1, by = 1/(K_IE+1))
    quantile_seq_IE = quantile_seq_IE[-c(1,K_IE+2)]
    
    ie_knots = quantile(X[,ind], quantile_seq_IE)
    
    IE_knots_stor[ind, ] = ie_knots
    
    ie_spl = bSpline(x = X[,ind], knots = ie_knots, intercept = FALSE)
    
    final_Xmat_IE = ie_spl
    
    IE_list[,,ind] = final_Xmat_IE
    
    S_IE_inv[,,ind] = penmatt(nspl_IE) 
    
  }
  
  ME_mat_MEs = NULL
  for(ind in 1:p)
  {
    
    ME_mat_MEs = cbind(ME_mat_MEs, ME_list[,,ind])
    
  }
  
  ME_mat = cbind(rep(1,n), Z, ME_mat_MEs)
  
  map_k_to_uv = matrix(0, nrow = choose(p,2), ncol = 3)
  
  for(k in 1:choose(p,2))
  {
    
    #Obtain inverse index (u,v)
    
    quad_dis = (2*p - 1)^2 - 8*k
    u = ceiling(0.5*((2*p-1) - quad_dis^0.5))
    v = p + k - (u*(p - 0.5*(u+1)))
    
    map_k_to_uv[k,1] = k
    map_k_to_uv[k,2] = u
    map_k_to_uv[k,3] = v
    
  }
  
  map_k_to_uv = map_k_to_uv - 1
  
  SigmaME = solve(penmatt(nspl_ME))
  SigmaME_inv = penmatt(nspl_ME)
  SigmaInt = solve(penmatt(nspl_IE))
  SigmaInt_inv = penmatt(nspl_IE)
  
  # SigmaME = diag(nspl_ME)
  # SigmaME_inv = diag(nspl_ME)
  # SigmaInt = diag(nspl_IE)
  # SigmaInt_inv = diag(nspl_IE)
  
  init_values = SIMinit(y = y,
                        X = X, 
                        Z = Z,
                        ME_list = ME_list, 
                        IE_list = IE_list, 
                        zero_ind = zero_ind, 
                        nugget = 0.01,
                        me_integral_constraint = me_integral_constraint)
  
  print(noquote(paste("########## Sampling initiated with MC = ", MC, " ########## ", sep = "")))
  
  SIM_model = SIDsampler_draws_adaptive_optimized(y, 
                                                  ME_mat, 
                                                  IE_list,
                                                  eps_MALA, 
                                                  c_HMC, 
                                                  L_HMC, 
                                                  MC,
                                                  n, 
                                                  p, 
                                                  p_cov,
                                                  SigmaME, 
                                                  SigmaME_inv,
                                                  SigmaInt, 
                                                  SigmaInt_inv,
                                                  nspl_ME, 
                                                  nspl_IE,
                                                  cutoff,
                                                  map_k_to_uv,
                                                  zero_ind,
                                                  accept_low,
                                                  accept_high,
                                                  accept_scale,
                                                  a_lamb,
                                                  b_lamb,
                                                  init_values,
                                                  precond)
  
  print(noquote(paste("########## Sampling completed with MC = ", MC, " ########## ", sep = "")))
  
  return(list("SIM_model" = SIM_model, 
              "ME_list" = ME_list,
              "IE_list" = IE_list,
              "ME_knots" = ME_knots_stor,
              "IE_knots" = IE_knots_stor,
              "data" = list("y" = y, "X" = X, "Z" = Z, "n" = n, "p" = p, "MC" = MC),
              "init_values" = init_values))
  
}
