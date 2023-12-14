sortedkMeans<-function(v)
{
  
  unsortedkMeans = kmeans(x = v, centers = 2, iter.max = 1)
  sortedkMeans = kmeans(x = v, centers = as.numeric(sort(unsortedkMeans$centers)), 
                        iter.max = 1)
  
  return(sortedkMeans)
  
}

S2MCluster<-function(v, cutoff)
{
  
  ## Sequentially 2-means cluster the elements of v into 2 clusters
  ## Stop whem M-m \leq cutoff, M = bigger mean
  ## Follow ``Variable selection using shrinkage priors'' by Li, Pati (2017)
  
  d = length(v)
  D = rep(0, d)
  
  ## Step 1: perform 2-means on v
  
  clust_int = sortedkMeans(v)
  clust1_ind = which(clust_int$cluster == 1) ##already in original scale
  
  D[clust1_ind] = 1
  D[-clust1_ind] = 2
  
  m = min(clust_int$centers)
  M = max(clust_int$centers)
  
  ## Sequential clustering
  
  #iter = 1
  
  while((M-m > cutoff) & (length(clust1_ind) > 2))
  {
    
    ## Update D to be indices in cluster with mean m
    
    D[clust1_ind] = 1
    D[-clust1_ind] = 2
    
    ## Perform 2-means only on the indices in smaller cluster
    
    clust_int = sortedkMeans(v[clust1_ind])
    
    ## Update m and M
    
    m = min(clust_int$centers)
    M = max(clust_int$centers)
    
    ## Translate to original indexing
    
    clust1_ind = clust1_ind[which(clust_int$cluster == 1)]
    
    ## Number of iterations
    
    #iter = iter + 1
    
  }
  
  return(D)
  
}

cutoffgridSignals<-function(v, cutoff_grid)
{
  
  L = length(cutoff_grid)
  numSignal = rep(0, L)
  
  for(l in 1:L)
  {
    
    cutoff = cutoff_grid[l]
    estD = S2MCluster(v, cutoff_grid[l])
    
    numSignal[l] = length(estD) - length(which(estD == 1))
    
  }
  
  return(numSignal)
  
}

cutoffDetect<-function(numSignal, cutoff_grid)
{
  
  unique_vals = unique(numSignal)
  L = length(unique_vals)
  
  unique_vals_diff = unique_vals[1:(L-1)] - unique_vals[2:L]
  
  max_ind = which.max(unique_vals_diff)
  cutoff_ind = min(which(numSignal == unique_vals[max_ind + 1]))
  
  return(cutoff_grid[cutoff_ind])
  
}

sparsifyVec<-function(v, length.grid = 1000, lower.thres = 0.00001)
{
  
  ## Induce hard zeroes in vector v \geq 0.
  ## Use S2M clustering idea in Li and Pati (2017).
  
  cutoff_grid = seq(min(v)+lower.thres, max(v), length.out = length.grid)
  numSignalVals = cutoffgridSignals(v = v, cutoff_grid = cutoff_grid)
  optimal_cutoff = cutoffDetect(numSignalVals, cutoff_grid)
  
  ## Run the model with chosen optimal cutoff.
  
  optimal_ind = S2MCluster(v, optimal_cutoff)
  
  ## Set indices for smaller cluster to zero.
  
  v[which(optimal_ind == 1)] = 0
  
  ## Return sparsified v.
  
  return(v)
  
}

## Given matrix of MCMC samples of int_max, int_min
## PIPOutput returns PIPs, PSPs, and PAPs of all interactions
## using the S2M approach to automatically choose cutoff.

PIPOutput<-function(int_max, int_min, length.grid = 1000, lower.thres = 0.00001)
{
  
  nMCMC = dim(int_max)[1]
  nInt = dim(int_max)[2]
  
  int_max_sparse = matrix(0, nrow = nMCMC, ncol = nInt)
  int_min_sparse = matrix(0, nrow = nMCMC, ncol = nInt)
  
  ## Equivalent to Single cutoff for all effects
  
  combined_int_sparse_vec = sparsifyVec(c(c(int_max), c(int_min)),
                                        length.grid = length.grid,
                                        lower.thres = lower.thres)
  dimVec = length(combined_int_sparse_vec)
  int_max_sparse_vec = combined_int_sparse_vec[1:(dimVec/2)]
  int_min_sparse_vec = combined_int_sparse_vec[((dimVec/2) + 1):dimVec]
  
  int_max_sparse = matrix(int_max_sparse_vec, nrow = nMCMC, ncol = nInt)
  int_min_sparse = matrix(int_min_sparse_vec, nrow = nMCMC, ncol = nInt)
  
  ## Evaluate posterior probabilities after sparsifying
  
  PIP_stor = rep(0, dim(int_max)[2])
  PSP_stor = rep(0, dim(int_max)[2])
  PAP_stor = rep(0, dim(int_max)[2])
  
  for(k in 1:dim(int_max)[2])
  {
    
    zero_int_ind = which((int_max_sparse[,k] == 0) & (int_min_sparse[,k] == 0))
    PIP_stor[k] = 1 - (length(zero_int_ind) / (dim(int_max)[1]))
    
    syn_int_ind = which((int_max_sparse[,k] > 0) & (int_min_sparse[,k] == 0))
    PSP_stor[k] = length(syn_int_ind) / (dim(int_max)[1])
    
    ant_int_ind = which((int_max_sparse[,k] == 0) & (int_min_sparse[,k] > 0))
    PAP_stor[k] = length(ant_int_ind) / (dim(int_max)[1])
    
  }
  
  outputList<-list("PIPInt" = PIP_stor,
                   "PSPInt" = PSP_stor,
                   "PAPInt" = PAP_stor)
  
  return(outputList)
  
}