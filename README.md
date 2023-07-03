# simulation2.R
(for data generation)
input parameters:
	B <- 200 # No. of simulation repetition
	n <- 400 # sample size
	pn <- 400 # No. of predictors
	prior_pi <- 0.4 # prior probabilities for two populations prior_pi and 1-prior_pi
	K <- 50 # No. of basis
	ck <- readRDS('E:/PhD/assistant works/functional classification/code/ck.rds') # keep fix for all simulations !!
	m <- 100 # spaced times
	phi <- 0 # measurement errors
	tau <- 0 # covariance structure difference between two populations
	rho <- 0.2 # controls the correlation among the functional predictors
	delta <- 2 # controls the signal strength
	rns <- c(1,3,5) # corresponding to the subsettings in each setting
	
	ifGaussian <- FALSE # FALSE - use non-Gaussian ksi_tilde
	df <- 3 # for non-Gaussian case

*	makeCluster(No. of multiprocess)

outputs:
	file_name_X # RDS data - dimension [400, 400, 100] - [n, pn, m] containing training and testing sets
	file_name_Y # RDS data - dimension [400]
	
# fclassification.R
(for classification of simulated data)
input parameters:
	B <- 200
	
	setting <- 8 # setting number(in total 8 general settings)
	delta <- 2 # same above
	rn <- 5 # same above
	phi <- 0 # same above
	tau <- 0 # same above
	(these 5 parameters are encoded in the file name of simulated data)
	
*	makeCluster(No. of multiprocess)

outputs:
	file_name_evatxt # evaluation table
	file_name_w_opt # optimal w
	file_name_w_select # optimal w before last projection step, used for calculating FNR FPR

# fclassification_realdata.R
(for classification of real data with 5-fold, 10-fold and 20-fold)
input parameters:
	file_name # pre-processed real data
	K_fold # No. of k-fold
	B # No. of repetition

*	makeCluster(No. of multiprocess)

outputs:
	(same as above)

# fclassification_realdata_leaveoneout.R
(for classification of real data with leave-one-out)
input parameters:
	file_name # pre-processed real data
	K_fold <- 121 # leave-one-out
	B # No. of repetition

*	makeCluster(No. of multiprocess)

outputs:
	(same as above)