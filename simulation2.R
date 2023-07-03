library(foreach)
library(doParallel)
start_time <- Sys.time()
B <- 200
n <- 400 # sample size
pn <- 400 # No. of predictors
prior_pi <- 0.4 # prior probabilities for two populations prior_pi and 1-prior_pi
K <- 50 # No. of basis
ck <- readRDS('E:/PhD/assistant works/functional classification/code/ck.rds') # keep fix for all simulations !!
m <- 100 # spaced times
phi <- 0 # measurement errors
tau <- 0 # covariance structure difference between two populations
rho <- 0.2 # controls the correlation among the functional predictors
delta <- 2
rns <- c(1,3,5)
ifGaussian <- FALSE
df <- 3 # for non-Gaussian case
for(rn in rns){

T_star <- c(1:rn)
k_lo <- 4 # for beta*
k_hi <- 6 # for beta*
non_zero_beta_coef <- c(k_lo:k_hi)

for (b in 1:B){
  Y <- rbinom(n, 1, 1-prior_pi) + 1
  
  beta_star <- matrix(0, rn, K) # parametrized vector in r^50*rn
  # calculate beta*
  for (j in T_star){
    beta_star[j, non_zero_beta_coef] <- sapply(non_zero_beta_coef, function(k){
      delta * (-1)**(k+1) * k**(2/5)
    })
  }
  beta_star_unstack <- rep(0, rn*K)
  for (j in 1:rn) {
    beta_star_unstack[((j-1)*K+1):(j*K)] <- beta_star[j,]
  }
  
  # orthonormal Fourier basis
  fBasis_eval <- matrix(1, K, m)
  fBasis_eval[1,] <- rep(1, m) 
  for (k in 2:K) {
    # create Fourier basis
    if (k %% 2 == 1){
      l <- (k + 1) / 2
      fBasis_eval[k,] <- sqrt(2) * sin((l - 1) * pi * (2*seq(0,1,1/(m-1)) - 1))
    }else{
      l <- k / 2
      fBasis_eval[k,] <- sqrt(2) * cos(l * pi * (2*seq(0,1,1/(m-1)) - 1))
    }
  }
  
  # simluate X(t)
  X <- array(0, dim = c(n, pn, m))
  
  simulateX <- function(i){
    y <- Y[i]
    epsi <- matrix(rnorm(pn*m), pn, m)
    
    
    # simulate ksi_tilde
    ksi_tilde <- matrix(0, nrow = pn, ncol = K)
    if(ifGaussian) {
      for (k in 1:K) {
        ksi_tilde[,k] <- rnorm(pn, mean = 0, sd = ck[k] / k**2)  # Gaussian scenario
      }
    }else {
      for (k in 1:K) {
        ksi_tilde[,k] <- sqrt(ck[k] / k**2 * (df-2) / df) * rt(pn, df)  # non-Gaussian scenario
      }
    }
    
    
    j2_bound <- rn # omit unnecessary calculation !!
    sigma_rho <- matrix(0, pn*K, j2_bound*K)
    # calculate sigma_rho if Y == 2
    rho_func_sum <- matrix(0, pn, j2_bound)
    if (y == 2){
      for (j1 in 1:pn){
        for (j2 in 1:j2_bound){
          rho_func <- function(j_prim) {
            rho ** (abs(j1 - j_prim) + abs(j2 - j_prim))
          }
          rho_func_sum[j1, j2] <- sum(sapply(c(1:pn), rho_func))
        }
      }
      
      var_coef <- diag(sapply(c(k_lo:k_hi), function(k){ck[k]/k**2}))
      for (j1 in 1:pn){
        for (j2 in 1: j2_bound){
          # s = k_hi-k_lo+1
          # sigma_rho[((j1-1)*K+k_lo):((j1-1)*K+k_hi), ((j2-1)*K+k_lo):((j2-1)*K+k_hi)] = matrix(1,s,s) * rho_func_sum[j1, j2] * var_coef
          sigma_rho[((j1-1)*K+k_lo):((j1-1)*K+k_hi), ((j2-1)*K+k_lo):((j2-1)*K+k_hi)] = rho_func_sum[j1, j2] * var_coef
        }
      }
    }
    
    # calculate ksi
    ksi <- matrix(0, pn, K)
    if (y == 2){
      for(j in 1:pn){
        ksi[j,] = sapply(c(1:K), function(k) {
          (1 + tau)**0.5 * ksi_tilde[,k] %*% rho**(abs(j - c(1:pn))) +
            sigma_rho[K*(j-1)+k,1:(j2_bound*K)] %*% beta_star_unstack
        })
      }
    } else {
      for(j in 1:pn){
        ksi[j,] = sapply(c(1:K), function(k) {
          ksi_tilde[,k] %*% rho**(abs(j - c(1:pn))) 
        })
      }
    }
    return(list('x' = ksi %*% fBasis_eval + phi * epsi))
    # cat('Xi(t)|Y, i=',i,' Y=',y,' complete.\n')
  }
  
  cl<- makeCluster(10)
  registerDoParallel(cl)
  x_tmp_list <- foreach(i = 1:n,
                        .combine = 'c') %dopar% simulateX(i)
  stopCluster(cl)
  for (i in 1:n){
    X[i,,] <- x_tmp_list[[i]]
  }
  
  # outputs
  cat('rn=',rn,',b=',b,'complete.\n')
  file_name_X <- paste0('G:/fclassification/sim_result_s8_d',as.character(delta),'/sim_tau',as.character(tau),'_phi',as.character(phi),'_beta2_nongaussian_delta',as.character(delta),'_rn',as.character(rn),'_x_b', as.character(b),'.rds')
  saveRDS(X,file = file_name_X)
  file_name_Y <- paste0('G:/fclassification/sim_result_s8_d',as.character(delta),'/sim_tau',as.character(tau),'_phi',as.character(phi),'_beta2_nongaussian_delta',as.character(delta),'_rn',as.character(rn),'_y_b', as.character(b),'.rds')
  saveRDS(Y,file = file_name_Y)
}# for B end

}# for rns end
end_time <- Sys.time()
print(end_time - start_time)
