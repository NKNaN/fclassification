library(stats)
library(fdapace)
library(glmnet)
library(foreach)
library(doParallel)
library(stringr)
library(reshape2)

set.seed(43)
y_df <- rep(0,122)
X_df <- array(0, dim = c(122, 64, 256))
file_name <- 'E:/PhD/assistant works/functional classification/code/real_data_preprocess/'
gz <- list.files(file_name)
gz <- sample(gz)
for (i in 1:length(gz)) {
  if (str_sub(gz[i], 4, 4) == 'c'){
    y_df[i] <- 1
  } else {
    y_df[i] <- 2
  }
  X_df[i,,] <- t(read.table(paste0(file_name,gz[i])))
}


K_fold <- 121
B <- 400
for(b in 1:B){

index <-  sort(sample(nrow(X_df), nrow(X_df)/K_fold))
X_test <- array(0,c(length(index),64,256))
for(i in length(index)){
  X_test[i,,] <- X_df[index[i],,]
}
y_test <- y_df[index]
X_train <-  X_df[-index,,]
y_train <- y_df[-index]
X <- X_train
y <- y_train

cat('b=',b)
# smooth

pooledcov <- function(Xlist, t){
  
  nx=lapply(Xlist,function(x){
    nrow(x)
  })
  n=sum(unlist(nx))
  
  outlist=vector("list",length(Xlist))
  for (j in 1:length(Xlist)){
    LX0 <- lapply(seq_len(ncol(t(Xlist[[j]]))), function(i) Xlist[[j]][i, ])
    Lt0 <- rep(list(t), nx[[j]])
    outlist[[j]]=FPCA(LX0,Lt0,optns  = list('useBinnedData' = 'OFF'))$fittedCov
  }
  
  
  Cov=0
  for (ii in 1:length(outlist)){
    Cov <- Cov+ (nx[[ii]] - 1) * outlist[[ii]]
    Cov <- Cov / (n - length(outlist) + 1)
  }
  
  
  return(list(Cov=Cov, singlecov=outlist))
}

getPhi <- function(Cov, FVEthreshold=99.99, t){
  eig <- eigen(Cov) #return a list with value and vector
  positiveInd <- eig[['values']] >= 0 #why not use dollar sign?
  if (sum(positiveInd) == 0) {
    stop('All eigenvalues are negative. The covariance estimate is incorrect.')
  }
  d <- eig[['values']][positiveInd] #same?
  
  eigenV <- eig[['vectors']][, positiveInd, drop=FALSE] #extract positive eigenvalue column
  
  
  FVE <- cumsum(d) / sum(d) * 100  # cumulative FVE for all available eigenvalues from fitted cov
  
  cutPoints <- min(which(FVE >= FVEthreshold)) # final number of component chosen based on FVE.Only 29 is choosen in demo.
  
  maxK <- cutPoints #maximum K candidate
  
  d <- d[1:maxK]
  
  proportion <- d/sum(d) #Contribution of each PC
  
  eigenV <- eigenV[, 1:maxK, drop=FALSE]
  
  
  # normalization
  muWork = 1:dim(eigenV)[1] #return the row number of eigenV. Since column need to be reduced
  
  
  phi0 <- apply(eigenV, 2, function(x) { #eigen V is matrix. 2 means apply function to column.
    x <- x / sqrt(fdapace::trapzRcpp(t, x^2))# divide each column by a constant. This integration see paper
    if ( 0 <= sum(x*muWork) )
      return(x)
    else
      return(-x)
  })
  
  return(list(phi = phi0, FVE = FVE, cutPoints = cutPoints, d = d, proportion = proportion))
}

getXi <- function(X.curve, phi, t){ #This function get the projection score. See the paper algorithm part.
  
  if(! is.numeric(t)) stop("t should be recorded time points!")
  
  
  xi <- NULL
  
  if(is.vector(X.curve)){
    xi <- apply(apply(phi, 2, function(x) {x * X.curve}), 2, function(x){return(trapzRcpp(t, x))})
  }
  
  else{
    for (i in 1:ncol(phi)){
      xi <- cbind(xi ,apply(t(matrix(rep(phi[, i], nrow(X.curve)), length(t))) * X.curve, 1, function(x){return(trapzRcpp(t, x))}))
    } #The Jth column is the series of Jth scores
  }
  return(xi)
}

# vector of Lkrj for each Xir(s)
Lkr <- function (xi, tau_list, tau_level, qlist, k) {
  q <- qlist[[as.character(k)]]
  tau <- tau_list[tau_level]
  quan <- q[tau_level,]
  res <- abs(xi - quan) * (mapply(function (x, y) {
    if (x > y) tau else (1 - tau)
  }, xi, quan))
  
  return(res)
}

z <- matrix(0, nrow(X), ncol(X))
z_test <- matrix(0, nrow(X_test), ncol(X_test))

funobj_train <- function(r) {
  Xr <- X[,r,]
  Xr_test <- array(0, dim=c(dim(X_test)[1],dim(X_test)[3]))
  Xr_test[1,] <- X_test[,r,]
  m <- ncol(Xr)
  t <- seq(0, 1, 1/(m-1))
  smoothlist <- apply(Xr, 1, function(x){  #use apply() to smooth
    spline(t,x,method = "natural")
  })
  smoothlist_test <- apply(Xr_test, 1, function(x){  #use apply() to smooth
    spline(t,x,method = "natural")
  })
  t <- smoothlist[[1]]$x #replace original t with smoothed t
  Xr <- matrix(0, nrow = nrow(Xr), ncol = length(t))
  Xr_test <- matrix(0, nrow = nrow(Xr_test), ncol = length(t))
  for (i in 1:nrow(Xr)){
    Xr[i, ] <- smoothlist[[i]]$y
  } #repalce original X with smoothed X
  for (i in 1:nrow(Xr_test)){
    Xr_test[i, ] <- smoothlist_test[[i]]$y
  }
  
  K <- length(unique(y))
  
  Xlist <- vector("list", K)
  bindX <- cbind(Xr, y)
  
  #construct index list to record original sequence of labels
  ind <- vector("list",K)
  for (j in 1:K){
    ind[[j]] <- which(bindX[, ncol(bindX)]==(unique(y))[j]) #record row index for each group i
  }
  names(ind) <- c(unique(y)[1: K])
  
  for (i in 1:K){
    Xlist[[i]] <- bindX[y==(unique(y))[i],]
    Xlist[[i]] <- matrix(as.numeric(Xlist[[i]][, -ncol(bindX)]),nrow(Xlist[[i]]))
  }
  names(Xlist) <- c(unique(y)[1:K])
  
  Covlist=pooledcov(Xlist,t)
  
  #Compute basis function by eigendecomposition
  FVEthreshold=99.99
  temp.phi <- getPhi(Covlist$Cov, FVEthreshold, t)
  
  fitphi <- temp.phi$phi
  
  # eigenV <- temp.phi$eigenV
  # 
  # cutPoints <- temp.phi$cutPoints
  # 
  # proportion <- temp.phi$proportion
  
  # compute xi
  xilist=vector("list",length(Xlist))
  for (jj in 1:length(Xlist)){
    xilist[[jj]] <- getXi(Xlist[[jj]], fitphi[, 1:temp.phi$cutPoints], t)
  }
  names(xilist) <- c(unique(y)[1:K])
  J=ncol(xilist[[1]])
  xilist_test <- getXi(Xr_test, fitphi[, 1:temp.phi$cutPoints], t)
  
  #compute quantile for each column(projection scores).
  qlist=vector("list",length(Xlist))
  a0=0.02
  M0=10
  theta=seq(a0, 1-a0, (1-2*a0)/(M0-1))
  Ltheta=(1:length(theta))
  
  for (kk in 1:length(Xlist)){ #quantile matrix for each group
    qlist[[kk]]=apply(xilist[[kk]], 2, function(x){ #replace all loops by apply()
      quantile(x,theta)
    })
  }
  names(qlist) <- c(unique(y)[1:K])
  
  # combine xi and ind
  xi_mat <- c() # total xi without grouping
  for (kk in 1:length(Xlist)){
    rownames(xilist[[kk]]) <- ind[[kk]]
    xi_mat <- rbind(xi_mat, xilist[[kk]])
  }
  
  # reorder xi_mat
  xi_mat <- sapply(1:nrow(Xr), function(i) {xi_mat[which(rownames(xi_mat) == i),]})
  xi_mat <- t(xi_mat)
  
  # compute Lkr and dir for xi
  dr <- rep(0, nrow(xi_mat))
  for (i in 1:nrow(xi_mat)){
    dr[i] <- sum(sapply(1:M0, function(m) {
      Lkr(xi_mat[i,], theta, m, qlist, 1) - Lkr(xi_mat[i,], theta, m, qlist, 2)}))
  }
  
  # compute Q_test for each dimension Xr_test
  dr_test <- rep(0, nrow(xilist_test))
  for (i in 1:nrow(xilist_test)){
    dr_test[i] <- sum(sapply(1:M0, function(m) {
      Lkr(xilist_test[i,], theta, m, qlist, 1) - Lkr(xilist_test[i,], theta, m, qlist, 2)}))
  }
  
  return(list('dr' = dr, 'dr_test' = dr_test))
}

z_tmp_list <- NULL
cl<- makeCluster(8)      
registerDoParallel(cl)
z_tmp_list <- foreach(r = 1:length(X[1,,1]),
        .combine = 'c',
        .packages = c('fdapace','stats')) %dopar% funobj_train(r)
stopCluster(cl)

for(i in 1:length(z_tmp_list)){
  if(i %% 2 == 1) {
    z[,(i+1)/2] <- z_tmp_list[[i]]
  } else {
    z_test[,i/2] <- z_tmp_list[[i]]
  }
}

# loss function
l2 <- function(y, w, z, lambd){
  n <- nrow(z)
  s <- sapply(1:n, function(i){
    if (y[i] == 1) log(pi_hat(w, z[i,])) else log(1 - pi_hat(w, z[i,]))
  })
  s <- sum(s) / n
  # l <- -s + lambd * sum(abs(w))
  return(s)
}

H <- 0.5 / min(sapply(1:nrow(z), function(i){norm(z[i,], type='2')}))
pi_hat <- function(w, zi) {
  res <- 1 / (1 + exp(-zi %*% w))
  return(res)
}
# p -> 1
# p -> 0
ci <- function(w, zi) {
  p <- max(min(pi_hat(w, zi), 0.0001), 0.9999)
  res <- p * (1 - p)
  return(res)
}
ri <- function(yi, w, zi) {
  p <- pi_hat(w, zi)
  zi %*% w + (yi - p) / ci(w, zi)
}
softhresg <- function(gamma, lambd) {
  sign(gamma) * max((abs(gamma) - lambd), 0)
}

projection <- function (y) {
  n <- length(y)
  y_order <- sort(y)
  for (i in (n-1):1) {
    
    ti <- (sum(y_order[(i+1):n]) - 1) / (n - i) # step 3
    
    if (ti >= y_order[i]) {
      return(sapply((y - ti), function(x) {max(x, 0)}))
    }
  }
  ti <- (sum(y) - 1) / n # step 4
  
  return(sapply((y - ti), function(x) {max(x, 0)}))
}

z<-z*H
z_test<-z_test*H

Y <- y - 1

projected_lasso <- function(z, y, tol_inner = 1e-6, tol_outer = 1e-4, lambd = 0.05, max_iter = 100){
  n <- nrow(z)
  w_outer <- rep(0, ncol(z))
  w0 <- rep(0, dim(z)[2])
  # fit <- glmnet(H*z, y, family = 'binomial', lambda = lambd/2)
  # w0 <- as.vector(fit$beta)
  
  w <- projection(w0) # w*
  proj_w <- projection(w0) # projected w*
  w_inner <- w # inner temp w record
  w_outer <- proj_w # outer temp storage of previous w0, ri ci computed from w_outer
  loss <- c()
  
  w_before <- NULL
  for(i in 1:max_iter){
    # inner loop - coordinate descent
    while(TRUE) {
      for(j in 1:ncol(z)) {
        ui <- sapply(1:nrow(z), function(i){ # ui computed from w which is updated in each coordinate descent
          ri(y[i], w_outer, z[i,]) - z[i,]%*%w
        })
        tmp <- sapply(1:nrow(z), function(i){
          ci(w_outer, z[i,]) * z[i,j] * ui[i]
        })
        w[j] <- softhresg(sum(tmp)/n, lambd)
      }
      
      w_diff_inner <- max(abs(w - w_inner))
      # cat('w_diff_inner: ', w_diff_inner,'\n')
      if (w_diff_inner < tol_inner) {
        break
      }
      w_inner <- w
    }
    
    w_before <- w
    # projection
    proj_w_new <- projection(w)
    # print(proj_w_new)
    loss <- c(loss, l2(y, proj_w_new, z, lambd))
    
    w_diff_outer <- max(abs(proj_w_new - proj_w))
    # cat('w_diff_outer: ', w_diff_outer, '\n')
    if (w_diff_outer < tol_outer) {
      break
    }
    w_outer <- proj_w_new
    proj_w <- proj_w_new
    w <- proj_w_new
    w_inner <- w
  }
  
  return(list('w_proj'=proj_w,'w_befporj'=w_before, 'loss'=loss))
}

# prediction to class '1' and '2'
predict <- function(d, w) {
  p <- length(w)
  Q <- d %*% w
  ret <- sapply(1:nrow(d), function(i){if (Q[i] > 0) 2 else 1}) # Q > 0, distance to group '1' larger than distance to group '2'
  return(ret)
}

evaluate <- function(y_true, y_pred) {
  mcr <- sum(y_true != y_pred)/length(y_true)
  tp <- 0
  fn <- 0
  fp <- 0
  tn <- 0
  for(i in 1:length(y_true)){
    if(y_true[i] == 1){
      if(y_pred[i] == 1){
        tp <- tp + 1
      }else{
        fn <- fn + 1
      }
    }else{
      if(y_pred[i] == 1){
        fp <- fp + 1
      }else{
        tn <- tn + 1
      }
    }
  }
  precision1 <- tp/(tp+fp)
  precision2 <- tn/(fn+tn)
  recall1 <- tp/(tp+fn)
  recall2 <- tn/(fp+tn)
  # if(is.na(precision1)){precision1=0}
  # if(is.na(precision2)){precision2=0}
  # if(is.na(recall1)){recall1=0}
  # if(is.na(recall2)){recall2=0}
  mp = (precision1 + precision2)/2
  mr = (recall1 + recall2)/2
  f1 = 2*mp*mr/(mp+mr)
  # if(is.na(f1)){f1=0}
  return(list('MCR' = mcr, 'Macro_P' = mp, 'Macro_R' = mr,
              'Macro_F1' = f1))
}

lambs = sapply(seq(0, 5, 1), function(x){0.0003*3**x})
nfolds=5
n <- nrow(z)
p <- ncol(z)
idx <- ceiling(sample(1:n)/n*nfolds) # nfolds intervals from 1/n to 1
mcr <- matrix(0, nfolds, length(lambs))

cv_inner <- function(i) {
  mcr <- rep(0, length(lambs))
  x_train <- z[idx!=i,]
  y_train <- Y[idx!=i]

  x_test  <- z[idx==i,]
  y_test  <- Y[idx==i]
  for (j in 1:length(lambs)) {
    fitted <- projected_lasso(x_train, y_train, lambd = lambs[j])
    y_pred <- predict(x_test, fitted$w_proj) - 1
    mcr[j] <- sum(y_test != y_pred)/length(y_test)
  }
  return(mcr)
}
cl<- makeCluster(8)
registerDoParallel(cl)
mcr <- foreach(i = 1:nfolds, .combine = 'rbind') %dopar% cv_inner(i)
stopCluster(cl)
cvm <- colMeans(mcr)
idx.min <- which.min(cvm)
lambd.min <- lambs[idx.min]

projected_lasso.fit <- projected_lasso(z, Y, lambd = lambd.min)
lasso_fit <- list('m' = mcr, 'cvm' = cvm, 'lambd.min' = lambd.min, 'fit' = projected_lasso.fit, 'lambd' = lambs)

w_opt <- lasso_fit$fit$w_proj
metrics <- evaluate(y_test, predict(z_test, w_opt))
write.table(metrics, file = paste0('E:/PhD/assistant works/functional classification/code/real_result_',as.character(K_fold),'f/evaluation.txt'),sep = '\t',
            append = T, col.names = F, row.names = F)
saveRDS(w_opt, file = paste0('E:/PhD/assistant works/functional classification/code/real_result_',as.character(K_fold),'f/w',as.character(b),'.rds'))
saveRDS(lasso_fit$fit$w_befporj, file = paste0('E:/PhD/assistant works/functional classification/code/real_result_',as.character(K_fold),'f/w_select_',as.character(b),'.rds'))
}