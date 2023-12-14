library(ks)

kde_cdf<-function(q, dat_vec, h)
{
  
  dat_vec = as.numeric(dat_vec)
  n = length(dat_vec)
  
  cdf_vec = pnorm(q = q, mean = dat_vec, sd = h)
  
  return(mean(cdf_vec))
  
}

trans_to_uniform<-function(X)
{
  
  X = as.matrix(X)
  n = dim(X)[1]
  p = dim(X)[2]
  
  new_X = matrix(0, nrow = n, ncol = p)
  all_bw = rep(0, p)
  
  for(j in 1:p)
  {
    
    dat_vec = X[,j]
    kde_mod = kde(x = dat_vec)
    
    kde_bw = kde_mod$h
    
    for(i in 1:n)
    {
      
      new_X[i,j] = kde_cdf(q = X[i,j], dat_vec = dat_vec, h = kde_bw)
      
    }
    
    all_bw[j] = kde_bw
    
    
  }
  
  res = list("new_X" = new_X, "all_bw" = all_bw)
  
  return(res)
  
}

map_k_to_uv<-function(k, p)
{
  
  quad_dis = (2*p - 1)^2 - 8*k
  u = ceiling(0.5*((2*p-1) - quad_dis^0.5))
  v = p + k - (u*(p - 0.5*(u+1)))
  
  return(c(u, v))
  
}

map_uv_to_k<-function(vec, p)
{
  
  u = vec[1]
  v = vec[2]
  
  k = u*(p - 0.5*(u+1)) + v - p
  
  return(k)
  
}

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

# SIM_extract<-function(whole_model_SIM, Xtest, good_index)
# {
#   
#   #### Returns whole surface and interaction estimate values ####
#   
#   SIM_model = whole_model_SIM$SIM_model
#   
#   MC = whole_model_SIM$data$MC
#   #good_index = (MC - 5000):MC
#   
#   ntest = dim(Xtest)[1]
#   p = dim(Xtest)[2]
#   
#   ME_knots = whole_model_SIM$ME_knots
#   IE_knots = whole_model_SIM$IE_knots
#   
#   nspl_ME = dim(ME_knots)[2] + 4
#   nspl_IE = dim(IE_knots)[2] + 3
#   
#   ## Construct basis function matrices for test data
#   
#   ME_test_mat = array(0, dim = c(ntest, nspl_ME, p))
#   IE_test_mat = array(0, dim = c(ntest, nspl_IE, p))
#   
#   for(j in 1:p)
#   {
#     
#     ME_knots_j = whole_model_SIM$ME_knots[j,]
#     IE_knots_j = whole_model_SIM$IE_knots[j,]
#     
#     ME_base_mat = bSpline(x = Xtest[,j], 
#                                knots = ME_knots_j,
#                                intercept = TRUE)
#     
#     ME_test_mat[,,j] = sweep(ME_base_mat, 2, colMeans(ME_base_mat))
#     
#     IE_test_mat[,,j] = bSpline(x = Xtest[,j], 
#                                knots = IE_knots_j,
#                                intercept = FALSE)
#     
#   }
#   
#   augME_mat_test = NULL
#   for(j in 1:p)
#   {
#     
#     augME_mat_test = cbind(augME_mat_test, ME_test_mat[,,j])
#     
#   }
#   
#   ## Extract whole surface estimates for test data and store error
#   
#   whole_sur_samples = matrix(0, nrow = ntest, ncol = length(good_index))
#   int_sur_samples = array(0, dim = c(ntest, length(good_index), choose(p,2)))
#   
#   int_storage = matrix(0, nrow = ntest, ncol = choose(p,2))
#   
#   for(m0 in 1:length(good_index))
#   {
#     
#     m = good_index[m0]
#     
#     intercept_sample = SIM_model$intercept[m]
#     ME_sample = augME_mat_test %*% SIM_model$ME_coeff[m,]
#     
#     for(k in 1:choose(p,2))
#     {
#       
#       u = map_k_to_uv(k, p)[1]
#       v = map_k_to_uv(k, p)[2]
#       
#       pos_IE_sample_k = ((IE_test_mat[,,u] %*% SIM_model$IE_pos_coeff_part1[m,,k])^2 * 
#                          (IE_test_mat[,,v] %*% SIM_model$IE_pos_coeff_part2[m,,k])^2)
#       neg_IE_sample_k = ((IE_test_mat[,,u] %*% SIM_model$IE_neg_coeff_part1[m,,k])^2 * 
#                          (IE_test_mat[,,v] %*% SIM_model$IE_neg_coeff_part2[m,,k])^2)
#       
#       IE_sample_k = pos_IE_sample_k - neg_IE_sample_k
#       
#       int_sur_samples[,m0,k] = IE_sample_k
#       
#       int_storage[,k] = IE_sample_k
#       
#     }
#     
#     whole_sur_samples[,m0] = intercept_sample + ME_sample + apply(int_storage, 1, sum)
#     
#   }
#   
#   whole_sur_est = rowMeans(whole_sur_samples)
#   #int_est = apply(int_sur_samples, c(1,3), mean)
#   int_est_samples = int_sur_samples
#   
#   return(list("whole_sur_est" = whole_sur_est,
#               "int_est_samples" = int_est_samples))
#   
# }

filled.contour2 <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
                                                                          length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
                             ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
                             levels = pretty(zlim, nlevels), nlevels = 20, color.palette = function(n) hcl.colors(n, 
                                                                                                                  "YlOrRd", rev = TRUE), col = color.palette(length(levels) - 
                                                                                                                                                               1), plot.title, plot.axes, key.title, key.axes, asp = NA, 
                             xaxs = "i", yaxs = "i", las = 1, axes = TRUE, frame.plot = axes, 
                             ...) 
{
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
  layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
  par(las = las)
  mar <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  par(mar = mar)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
              yaxs = "i")
  rect(0, levels[-length(levels)], 1, levels[-1L], col = col,
       border = NA)
  if (missing(key.axes)) {
    if (axes) 
      axis(4)
  }
  else key.axes
  box()
  if (!missing(key.title)) 
    key.title
  mar <- mar.orig
  mar[4L] <- 1
  par(mar = mar)
  plot.new()
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  .filled.contour(x, y, z, levels, col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}