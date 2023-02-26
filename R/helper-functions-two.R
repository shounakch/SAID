library(bkmr)
# library(SynInt)
# source("Competitors_fcts.R")
# source("~/Research/Shape Constrained Inference/final-simulations/MixSelect-master/MixSelect.R", echo=TRUE)
# source("~/Research/Shape Constrained Inference/final-simulations/MixSelect-master/predict_MixSelect.R", echo=TRUE)

#Maps an integer k to (u,v) such that 1 \leq u < v \leq p.
map_k_to_uv<-function(k, p)
{
  
  quad_dis = (2*p - 1)^2 - 8*k
  u = ceiling(0.5*((2*p-1) - quad_dis^0.5))
  v = p + k - (u*(p - 0.5*(u+1)))
  
  return(c(u, v))
  
}

#Reverse mapping.
map_uv_to_k<-function(vec, p)
{
  
  u = vec[1]
  v = vec[2]
  
  k = u*(p - 0.5*(u+1)) + v - p
  
  return(k)
  
}

#Decomposes a vector into v = v_{+} - v_{-}.
#Also returns mean(v_{+}) and mean(v_{-}).
vec_decomp<-function(v)
{
  
  v1 = sapply(v, max, 0)
  v2 = sapply(-v, max, 0)
  
  int1 = mean(v1)
  int2 = mean(v2)
  
  mylist = list("vplus" = v1,
                "vminus" = v2,
                "intplus" = int1,
                "intminus" = int2)
  
  return(mylist)
  
}

##Extracts the intplus and intminus parts from a list of above vec_decomps.

intplus_extract<-function(L)
{
  
  output = L$intplus
  
  return(output)
  
}

intminus_extract<-function(L)
{
  
  output = L$intminus
  
  return(output)
  
}

vec_max<-function(v)
{
  
  return(sapply(v, max, 0))
  
}

############## SIM Evaluation ################

#Extracts whole surface and IE estimate posterior samples for SIM.
#Post-processing SIM function - add this to SynInt package.

SIM_extract<-function(whole_model_SIM, Xtest, good_index, me_integral_constraint)
{
  
  #### Returns whole surface and interaction estimate values ####
  ## Assumes there are no covariates. ##
  
  ###First construct bSpline models with which to predict###
  
  SIM_model = whole_model_SIM$SIM_model
  n = whole_model_SIM$data$n
  MC = whole_model_SIM$data$MC
  X = whole_model_SIM$data$X
  
  ntest = dim(Xtest)[1]
  p = dim(Xtest)[2]
  
  ME_knots = whole_model_SIM$ME_knots
  IE_knots = whole_model_SIM$IE_knots
  
  if(me_integral_constraint == TRUE)
  {
    
    nspl_ME = dim(ME_knots)[2] + 4
    nspl_IE = dim(IE_knots)[2] + 3
    
  }else
  {
    
    nspl_ME = dim(ME_knots)[2] + 3
    nspl_IE = dim(IE_knots)[2] + 3
    
  }
  
  #Construct ME and IE spline objects for extrapolation
  
  ME_list = array(0, dim = c(n, nspl_ME, p))
  IE_list = array(0, dim = c(n, nspl_IE, p))
  
  ME_int_adjust = matrix(0, nrow = nspl_ME, ncol = p)
  
  ind = 1
  ME_spline_obj_list = vector(mode = "list", length = p)
  IE_spline_obj_list = vector(mode = "list", length = p)
  
  Xwhole = rbind(X, Xtest)
  
  ## Construct basis function matrices for test data
  
  ME_test_mat = array(0, dim = c(ntest, nspl_ME, p))
  IE_test_mat = array(0, dim = c(ntest, nspl_IE, p))
  
  for(ind in 1:p)
  {
    
    ### For main effects ###
    
    if(me_integral_constraint == TRUE)
    {
      
      me_spl = bSpline(x = Xwhole[,ind],
                       knots = ME_knots[ind,],
                       intercept = TRUE)
      
      final_Xmat_ME = as.matrix(me_spl)
      final_Xmat_ME = sweep(final_Xmat_ME, 2, colMeans(final_Xmat_ME))
      
      ME_list[,,ind] = final_Xmat_ME[1:n,]
      ME_test_mat[,,ind] = final_Xmat_ME[-c(1:n),]
      
    }else
    {
      
      me_spl = bSpline(x = Xwhole[,ind],
                       knots = ME_knots[ind,],
                       intercept = FALSE)
      
      final_Xmat_ME = as.matrix(me_spl)
      
      ME_list[,,ind] = final_Xmat_ME[1:n,]
      ME_test_mat[,,ind] = final_Xmat_ME[-c(1:n),]
      
    }
    
    ME_spline_obj_list[[ind]] = me_spl
    
    ### For interaction effects ###
    
    ie_spl = bSpline(x = Xwhole[,ind],
                     knots = IE_knots[ind,],
                     intercept = FALSE)
    
    final_Xmat_IE = as.matrix(ie_spl)
    
    IE_list[,,ind] = final_Xmat_IE[1:n,]
    IE_test_mat[,,ind] = final_Xmat_IE[-c(1:n),]
    
    IE_spline_obj_list[[ind]] = ie_spl
    
  }
  
  # j = 1
  # 
  # for(j in 1:p)
  # {
  #   
  #   #ME_test_mat[,,j] = predict(ME_spline_obj_list[[j]], Xtest[,j])
  #   
  #   ME_test_mat[,,j] = as.matrix(ME_spline_obj_list[[j]])[-c(1:n),]
  #   
  #   #IE_test_mat[,,j] = predict(IE_spline_obj_list[[j]], Xtest[,j])
  #   
  #   IE_test_mat[,,j] = as.matrix(IE_spline_obj_list[[j]])[-c(1:n),]
  #   
  # }
  
  augME_mat_test = NULL
  for(j in 1:p)
  {
    
    augME_mat_test = cbind(augME_mat_test, ME_test_mat[,,j])
    
  }
  
  ## Extract whole surface estimates for test data and store error
  
  whole_sur_samples = matrix(0, nrow = ntest, ncol = length(good_index))
  int_sur_samples = array(0, dim = c(ntest, length(good_index), choose(p,2)))
  int_sur_train_samples = array(0, dim = c(n, length(good_index), choose(p,2)))
  
  for(k in 1:choose(p,2))
  {
    
    u = map_k_to_uv(k, p)[1]
    v = map_k_to_uv(k, p)[2]
    
    Theta1_k = SIM_model$IE_pos_coeff_part1[good_index,,k]
    Phi1_k = SIM_model$IE_pos_coeff_part2[good_index,,k]
    Theta2_k = SIM_model$IE_neg_coeff_part1[good_index,,k]
    Phi2_k = SIM_model$IE_neg_coeff_part2[good_index,,k]
    
    # Store at test points
    
    P1_samples_k = (IE_test_mat[,,u] %*% t(Theta1_k))^2
    P2_samples_k = (IE_test_mat[,,v] %*% t(Phi1_k))^2
    N1_samples_k = (IE_test_mat[,,u] %*% t(Theta2_k))^2
    N2_samples_k = (IE_test_mat[,,v] %*% t(Phi2_k))^2
    
    int_samples_k = (P1_samples_k * P2_samples_k) - (N1_samples_k * N2_samples_k)
    
    int_sur_samples[,,k] = int_samples_k
    
    # Store interaction samples at training points
    
    P1_samples_train_k = (IE_list[,,u] %*% t(Theta1_k))^2
    P2_samples_train_k = (IE_list[,,v] %*% t(Phi1_k))^2
    N1_samples_train_k = (IE_list[,,u] %*% t(Theta2_k))^2
    N2_samples_train_k = (IE_list[,,v] %*% t(Phi2_k))^2
    
    int_sur_train_samples[,,k] = (P1_samples_train_k * P2_samples_train_k) - 
      (N1_samples_train_k * N2_samples_train_k)
    
  }
  
  int_sum_samples = apply(int_sur_samples, c(1,2), sum)
  
  whole_sur_samples = (rep(1, dim(Xtest)[1]) %*% t(SIM_model$intercept[good_index])) +
                      (augME_mat_test %*% t(SIM_model$ME_coeff[good_index,])) +
                      int_sum_samples
  
  return(list("whole_sur_samples" = whole_sur_samples,
              "int_sur_samples" = int_sur_samples,
              "ME_spline_list" = ME_spline_obj_list,
              "IE_spline_list" = IE_spline_obj_list,
              "int_sur_train_samples" = int_sur_train_samples))
  
}

############### Non SIM-evaluation #######################

#Constructs a big matrix to obtain whole surface (H)
#estimates evaluated on test data for non-SIM methods.
whole_and_int_test<-function(Xtest)
{

  nt = dim(Xtest)[1]
  p = dim(Xtest)[2]
  K = choose(p,2)

  bigXtest1 = rep(0, p)

  bigXtest2 = Xtest

  bigXtest3 = NULL

  for(j in 1:p)
  {

    M_j = matrix(0, nrow = nt, ncol = p)
    M_j[,j] = Xtest[,j]
    bigXtest3 = rbind(bigXtest3, M_j)

  }

  bigXtest4 = NULL

  for(k in 1:K)
  {

    k_to_uv = map_k_to_uv(k, p)
    u = k_to_uv[1]
    v = k_to_uv[2]

    I_k = matrix(0, nrow = nt, ncol = p)

    I_k[,u] = Xtest[,u]
    I_k[,v] = Xtest[,v]

    bigXtest4 = rbind(bigXtest4, I_k)

  }

  bigXtest = rbind(bigXtest1, bigXtest2, bigXtest3, bigXtest4)

  return(bigXtest)

}

#Provides different components of global surface H estimated.
#intercept, ME, IEs
isolatefn<-function(bigvec, nt, p)
{

  K = choose(p,2)

  #Follow notation in SIM paper
  #First element is H(0_p)
  #Next nt elements are H evaluated at test points
  #Next nt * p elements are H evaluated at x_j * e_j, 1 \leq j \leq p
  #Next nt * K elements are H evaluated at x_u * e_u + x_v * e_v, 1 \leq u < v \leq p

  if(length(bigvec) == 1 + nt + (p*nt) + (K*nt))
  {

    H0 = 0
    Htest = rep(0, nt)
    Hmarginal = matrix(0, nrow = nt, ncol = p)
    Hjoint = matrix(0, nrow = nt, ncol = K)

    H0 = bigvec[1]

    Htest = bigvec[2:(nt + 1)]

    Hallmarginals = bigvec[(nt+2):(1 + nt + (p*nt))]
    Hmarginal = matrix(Hallmarginals, nrow = nt, ncol = p)

    Halljoints = bigvec[-c(1:(1 + nt + (p*nt)))]
    Hjoint = matrix(Halljoints, nrow = nt, ncol = K)

    output = list("H0" = H0,
                  "Htest" = Htest,
                  "Hmarginal" = Hmarginal,
                  "Hjoint" = Hjoint,
                  "nt" = nt,
                  "p" = p)

    return(output)

  }else
  {

    print("Error! Length of input incorrect!")

  }

}

#Extract interaction h(x_u, x_v) = H(x_u e_u + x_v e_v) - H(x_u e_u) -
#                                  H(x_v e_v) + H(0_{p}).
int_extract<-function(isolatevec)
{

  H0 = isolatevec$H0
  Hmarginal = isolatevec$Hmarginal
  Hjoint = isolatevec$Hjoint
  nt = isolatevec$nt
  p = isolatevec$p
  K = choose(p,2)

  int_stor = matrix(0, nrow = nt, ncol = K)

  for(k in 1:K)
  {

    k_to_uv = map_k_to_uv(k, p)
    u = k_to_uv[1]
    v = k_to_uv[2]

    int_stor[,k] = Hjoint[,k] - Hmarginal[,u] - Hmarginal[,v] + rep(H0, nt)

  }

  return(int_stor)

}

############## Frequentist Estimates ##########

# int_methods_freq<-function(name, y, X, Xtest)
# {
# 
#   # Compute estimates for interaction #k using 'name' method
#   # name \in {'Hiernet', 'PIE', 'FAMILY', 'RAMP'}
#   # kth interaction extracted. For p=2, k=1 ALWAYS.
# 
#   p = dim(X)[2]
#   nt = dim(Xtest)[1]
# 
#   fnc_name = paste(name, "_fct", sep = "")
# 
#   selected_fnc <- match.fun(fnc_name)
# 
#   # Construct the giant Xtest_new matrix to isolate interactions
# 
#   Xtest_new = whole_and_int_test(Xtest)
# 
#   fnc_mod = selected_fnc(y = y,
#                          X = X,
#                          X_test = Xtest_new)
# 
#   #### Get the interaction ####
# 
#   Heval = fnc_mod$y_pred
#   Hisolate = isolatefn(Heval, nt, p)
# 
#   whole_surface_est = Hisolate$Htest
#   est_int = int_extract(Hisolate)
# 
#   return(list("whole_surface_est" = whole_surface_est,
#               "interaction_est" = est_int))
# 
# }

######################## BKMR #########################

# int_methods_BKMR <- function(y, X, Xtest, iter)
# {
#   
#   BKMRmod = BKMR_fct(y = y,
#                      X = X,
#                      X_test = Xtest, 
#                      Z = NULL,
#                      iter = iter)
#   
#   whole_surface_est = BKMRmod$y_pred
#   
#   return(list("whole_surface_est" = whole_surface_est))
#   
# }
# 
# ################ MixSelect ################
# 
# int_methods_MixSelect <- function(y, X, Xtest,
#                                   MC = 5000)
# {
# 
#   MixSelectmod = MixSelect(y = y,
#                            X = X,
#                            nrun = MC,
#                            burn = MC - 1000)
# 
#   final_Xtest = whole_and_int_test(Xtest)
# 
#   MixSelectest = predict_MixSelect(MixSelectmod,
#                                    X_test = final_Xtest)
#   Heval = MixSelectest
# 
#   # Extract whole and interaction surfaces from Heval #
# 
#   nt = dim(Xtest)[1]
#   p = dim(Xtest)[2]
# 
#   Hisolate = isolatefn(Heval, nt, p)
# 
#   whole_surface_est = Hisolate$Htest
#   est_int = int_extract(Hisolate)
# 
#   return(list("whole_surface_est" = whole_surface_est,
#               "interaction_est" = est_int))
# 
# }

############## 

########################################################################

################### OLD CODE ###########################

# ##############################################################
# 

# 
# ################ BKMR ################
# 
# int_methods_BKMR <- function(y, X, Xtest,
#                              MC = 5000)
# {
#   
#   BKMRmod = bkmr::kmbayes(y = y,
#                           Z = X,
#                           iter = MC,
#                           varsel = TRUE)
#   
#   final_Xtest = whole_and_int_test(Xtest)
#   
#   BKMRest = ComputePostmeanHnew(BKMRmod, 
#                                 Znew = final_Xtest, 
#                                 method = "approx")
#   Heval = BKMRest$postmean
#   
#   # Extract whole and interaction surfaces from Heval #
#   
#   nt = dim(Xtest)[1]
#   p = dim(Xtest)[2]
#   
#   Hisolate = isolatefn(Heval, nt, p)
#   
#   whole_surface_est = Hisolate$Htest
#   est_int = int_extract(Hisolate)
#   
#   return(list("whole_surface_est" = whole_surface_est,
#               "interaction_est" = est_int))
#   
# }
# 
# ################ MixSelect ################
# 
# int_methods_MixSelect <- function(y, X, Xtest,
#                                   MC = 5000)
# {
#   
#   MixSelectmod = MixSelect(y = y,
#                            X = X,
#                            nrun = MC,
#                            burn = MC - 1000)
#   
#   final_Xtest = whole_and_int_test(Xtest)
#   
#   MixSelectest = predict_MixSelect(MixSelectmod,
#                                    X_test = final_Xtest)
#   Heval = MixSelectest
#   
#   # Extract whole and interaction surfaces from Heval #
#   
#   nt = dim(Xtest)[1]
#   p = dim(Xtest)[2]
#   
#   Hisolate = isolatefn(Heval, nt, p)
#   
#   whole_surface_est = Hisolate$Htest
#   est_int = int_extract(Hisolate)
#   
#   return(list("whole_surface_est" = whole_surface_est,
#               "interaction_est" = est_int))
#   
# }

# ######################## Frequentist Methods ##########################
# 
# int_methods_freq<-function(name, y, X, Xtest)
# {
# 
#   # Compute estimates for whole and interaction surfaces using 'name' method
#   # name \in {'Hiernet', 'PIE', 'FAMILY', 'RAMP'}
# 
#   p = dim(X)[2]
#   nt = dim(Xtest)[1]
# 
#   fnc_name = paste(name, "_fct", sep = "")
# 
#   selected_fnc <- match.fun(fnc_name)
# 
#   # Construct the giant Xtest_new matrix to isolate interactions
# 
#   fnc_mod = selected_fnc(y = y,
#                          X = X,
#                          X_test = Xtest)
#   
#   #### Get the whole surface estimate ####
#   
#   whole_surface_est = fnc_mod$y_pred
# 
#   #### Get the interaction ####
#   
#   int_surface_est = matrix(0, nrow = nt, ncol = choose(p,2))
#   int_coeffs = fnc_mod$Omega
#   
#   for(k in 1:choose(p,2))
#   {
#     
#     u = map_k_to_uv(k, p)[1]
#     v = map_k_to_uv(k, p)[2]
#     
#     int_surface_est[,k] = int_coeffs[u,v] * Xtest[,u] * Xtest[,v]
#     
#   }
#   
#   return(list("whole_surface_est" = whole_surface_est,
#               "interaction_est" = int_surface_est))
# 
# }








