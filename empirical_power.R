# Calculate empirical power for stepped wedge and staircase designs
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

library(MASS)
library(lme4)
library(glmmTMB)
library(pbkrtest)

source('functions_releff.R')

emp_power <- function(nsim, SWdes, K, m, ICC, CAC, theta){
  # Calculates empirical power based on nsim simulated trial datasets
  
  # Generate trial dataset and fit models, nsim times
  resrep <- replicate(nsim, fitmodels(SWdes, K, m, ICC, CAC, theta), simplify=FALSE)
  res <- list_rbind(resrep) # Convert from list to dataframe
  
  # Calculate rejection proportion for empirical power
  nsimSW <- sum(!is.na(res$estSW) & !is.na(res$seSW))
  pwrSW <- sum(abs(res$estSW)/res$seSW > 1.96, na.rm=TRUE)/nsimSW
  nsimSC <- sum(!is.na(res$estSC) & !is.na(res$seSC))
  pwrSC <- sum(abs(res$estSC)/res$seSC > 1.96, na.rm=TRUE)/nsimSC
  
  empvarSW <- var(res$estSW, na.rm=TRUE)
  empvarSC <- var(res$estSC, na.rm=TRUE)
  
  powvals <- data.frame(power_SW=pwrSW, power_SC=pwrSC,
                        emp_var_SW=empvarSW, emp_var_SC=empvarSC,
                        nsimSW=nsimSW, nsimSC=nsimSC)
  
  return(powvals)
}

fitmodels <- function(SWdes, K, m, ICC, CAC, theta){
  # Generates a single simulated trial dataset, fits corresponding model and
  # outputs the treatment effect estimate and standard error
  
  # Generate trial datasets
  dats <- gentrialdata(SWdes, K, m, ICC, CAC, theta)
  datSW <- dats[[1]]
  datSC <- dats[[2]]
  
  # Fit model
  if(CAC == 1){
    fitSW <- lmer(Y ~ treat + time + (1|cluster), data=datSW, REML=TRUE)
    estSW <- fixef(fitSW)['treat']
    seSW <- sqrt(vcov(fitSW)['treat','treat'])
    
    fitSC <- lmer(Y ~ treat + time + (1|cluster), data=datSC, REML=TRUE)
    estSC <- fixef(fitSC)['treat']
    seSC <- sqrt(vcov(fitSC)['treat','treat'])
  }else{
    fitSW <- fitBEmodel(datSW)
    estSW <- ifelse(is.null(fitSW), NA, fixef(fitSW)['treat'])
    seSW <- ifelse(is.null(fitSW), NA, sqrt(vcov(fitSW)['treat','treat']))
    
    fitSC <- fitBEmodel(datSC)
    estSC <- ifelse(is.null(fitSC), NA, fixef(fitSC)['treat'])
    seSC <- ifelse(is.null(fitSC), NA, sqrt(vcov(fitSC)['treat','treat']))
  }
  
  res1 <- data.frame(estSW=estSW, seSW=seSW, estSC=estSC, seSC=seSC)
  return(res1)
}


fitBEmodel <- function(dat){
  tryCatch(
    expr = {
      fit <- lmer(Y ~ treat + time + (1|cluster) + (1|clustper), data=dat, REML=TRUE)
    },
    error = function(e){
      message('Caught an error!')
      message(e)
      
      return(NULL)
    }
  )    
}


gentrialdata <- function(SWdes, K, m, ICC, CAC, theta){
  # Generates a single simulated trial dataset from a trial with given design
  # and trial configuration, for a discrete-time decay within-cluster
  # correlation structure with given correlation parameters and treatment effect
  #
  # Inputs:
  #   SWdes - stepped wedge design schematic with unique sequences only
  #   K - number of clusters assigned to each sequence in SWdes
  #   m - cluster-period size
  #   ICC - within-period intracluster correlation
  #   CAC - cluster autocorrelation (CAC=1 returns exchangeable correlation)
  # Output:
  #   One trial dataset, as a vector with length equal to the total number of
  #   measurements
  
  # Determine covariance parameter values based on given correlation parameters
  # Assuming total variance of 1: sig2T = 1 = sig2C + sig2e
  sig2C <- ICC
  sig2e <- 1 - sig2C
  
  # Determine remaining trial configuration parameters based on given design
  # matrix
  S <- dim(SWdes)[1] # number of unique sequences
  Tp <- dim(SWdes)[2] # total number of time periods (even if not all observed
  # in each sequence)
  
  fullSWdes <- SWdesmat(S, K)
  nclust <- dim(fullSWdes)[1] # number of clusters
  
  # Get vector of treatment effect indicators from given design matrix
  Xvec <- rep(as.vector(t(fullSWdes)), each=m)
  
  # Create vector of treatment effect indicators for embedded basic staircase
  # design matrix
  SCdes <- SCdesmat(S, K, 1, 1)
  XvecSC <- rep(as.vector(t(SCdes)), each=m)
  
  # Set up indices for cluster, cluster-period, and period
  clusterind <- rep(1:nclust, each=Tp*m)
  clustperind <- rep(1:(nclust*Tp), each=m)
  perind <- rep(1:Tp, each=m, times=nclust)
  
  # Time period effects (categorical)
  betavec <- perind
  
  # Random samples of each term
  # Block-exchangeable correlation
  exponent_matrix <- matrix(1, nrow=Tp, ncol=Tp) - diag(1, nrow=Tp, ncol=Tp)
  CP <- mvrnorm(n=nclust, mu=rep(0,Tp), Sigma=sig2C * CAC^(exponent_matrix))
  CPvec <- rep(as.vector(t(CP)), each=m)
  e <- rnorm(nclust*Tp*m, mean=0, sd=sqrt(sig2e))
  
  # Combine terms to get outcome data for SW design
  Y <- betavec + Xvec*theta + CPvec + e
  
  # Create data frame with everything needed for fitting the model
  datSW <- data.frame(Y=Y, cluster=as.factor(clusterind), time=as.factor(perind), clustper=as.factor(clustperind), treat=Xvec)
  
  # Get subset of data corresponding to embedded basic staircase design
  datSC <- datSW[!is.na(XvecSC),]
  
  dats <- list(datSW, datSC)
  return(dats)
}

set.seed(3004)

# Calculate empirical power for the designs under block-exchangeable correlation
res20_BE <- emp_power(1000, SWdesmat(5,1), 8, 20, 0.032, 0.93, 0.15)
res20_BE_typeI <- emp_power(1000, SWdesmat(5,1), 8, 20, 0.032, 0.93, 0)

res30_BE <- emp_power(1000, SWdesmat(5,1), 8, 30, 0.032, 0.93, 0.15)
res30_BE_typeI <- emp_power(1000, SWdesmat(5,1), 8, 30, 0.032, 0.93, 0)

res20_BE_50clust <- emp_power(1000, SWdesmat(5,1), 10, 20, 0.032, 0.93, 0.15)
res20_BE_50clust_typeI <- emp_power(1000, SWdesmat(5,1), 10, 20, 0.032, 0.93, 0)

res20_BE_60clust <- emp_power(1000, SWdesmat(5,1), 12, 20, 0.032, 0.93, 0.15)
res20_BE_60clust_typeI <- emp_power(1000, SWdesmat(5,1), 12, 20, 0.032, 0.93, 0)

# Empirical relative efficiency
res20_BE$emp_var_SW/res20_BE$emp_var_SC
res20_BE$emp_var_SW/res30_BE$emp_var_SC
res20_BE$emp_var_SW/res20_BE_50clust$emp_var_SC
res20_BE$emp_var_SW/res20_BE_60clust$emp_var_SC

