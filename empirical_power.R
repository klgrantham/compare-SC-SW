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
  
  # Type I error rate
  nsimSW0 <- sum(!is.na(res$estSW0) & !is.na(res$seSW0))
  typeISW <- sum(abs(res$estSW0)/res$seSW0 > 1.96, na.rm=TRUE)/nsimSW0
  nsimSC0 <- sum(!is.na(res$estSC0) & !is.na(res$seSC0))
  typeISC <- sum(abs(res$estSC0)/res$seSC0 > 1.96, na.rm=TRUE)/nsimSC0
  
  powvals <- data.frame(power_SW=pwrSW, power_SC=pwrSC,
                        emp_var_SW=empvarSW, emp_var_SC=empvarSC,
                        nsimSW=nsimSW, nsimSC=nsimSC,
                        typeI_SW=typeISW, typeI_SC=typeISC,
                        nsimSW0=nsimSW0, nsimSC0=nsimSC0)
  
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
    
    fitSW0 <- lmer(Y0 ~ treat + time + (1|cluster), data=datSW, REML=TRUE)
    estSW0 <- fixef(fitSW0)['treat']
    seSW0 <- sqrt(vcov(fitSW0)['treat','treat'])
    
    fitSC <- lmer(Y ~ treat + time + (1|cluster), data=datSC, REML=TRUE)
    estSC <- fixef(fitSC)['treat']
    seSC <- sqrt(vcov(fitSC)['treat','treat'])
    
    fitSC0 <- lmer(Y0 ~ treat + time + (1|cluster), data=datSC, REML=TRUE)
    estSC0 <- fixef(fitSC0)['treat']
    seSC0 <- sqrt(vcov(fitSC0)['treat','treat'])
  }else{
    fitSW <- fitDTDmodel(datSW)
    estSW <- ifelse(is.null(fitSW), NA, fixef(fitSW)[[1]]['treat'])
    seSW <- ifelse(is.null(fitSW), NA, getSE(fitSW))

    fitSW0 <- fitDTDmodel0(datSW)
    estSW0 <- ifelse(is.null(fitSW0), NA, fixef(fitSW0)[[1]]['treat'])
    seSW0 <- ifelse(is.null(fitSW0), NA, getSE(fitSW0))

    fitSC <- fitBEmodelSC(datSC)
    estSC <- ifelse(is.null(fitSC), NA, fixef(fitSC)['treat'])
    seSC <- ifelse(is.null(fitSC), NA, sqrt(vcov(fitSC)['treat','treat']))

    fitSC0 <- fitBEmodelSC0(datSC)
    estSC0 <- ifelse(is.null(fitSC0), NA, fixef(fitSC0)['treat'])
    seSC0 <- ifelse(is.null(fitSC0), NA, sqrt(vcov(fitSC0)['treat','treat']))
  }
  
  res1 <- data.frame(estSW=estSW, seSW=seSW, estSC=estSC, seSC=seSC,
                     estSW0=estSW0, seSW0=seSW0, estSC0=estSC0, seSC0=seSC0)
  return(res1)
}

fitDTDmodel <- function(dat){
  tryCatch(
    expr = {
      fitSW <- glmmTMB(Y ~ treat + time + ar1(time + 0|cluster), data=dat, REML=TRUE)
    },
    error = function(e){
      message('Caught an error!')
      message(e)
      
      return(NULL)
    }
  )
}

fitDTDmodel0 <- function(dat){
  tryCatch(
    expr = {
      fitSW0 <- glmmTMB(Y0 ~ treat + time + ar1(time + 0|cluster), data=dat, REML=TRUE)
    },
    error = function(e){
      message('Caught an error!')
      message(e)
      
      return(NULL)
    }
  )    
}

getSE <- function(fit){
  tryCatch(
    expr = {
      sqrt(vcov(fit)[[1]]['treat','treat'])
    },
    error = function(e){
      message('Caught an error!')
      message(e)
      
      return(NA)
    }
  )
}

fitBEmodelSC <- function(datSC){
  tryCatch(
    expr = {
      fitSC <- lmer(Y ~ treat + time + (1|cluster) + (1|clustper), data=datSC, REML=TRUE)
    },
    error = function(e){
      message('Caught an error!')
      message(e)
      
      return(NULL)
    }
  )    
}

fitBEmodelSC0 <- function(datSC){
  tryCatch(
    expr = {
      fitSC0 <- lmer(Y0 ~ treat + time + (1|cluster) + (1|clustper), data=datSC, REML=TRUE)
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
  # Discrete-time decay correlation
  CP <- mvrnorm(n=nclust, mu=rep(0,Tp),
                Sigma=sig2C * CAC^(as.matrix(dist(1:Tp))))
  CPvec <- rep(as.vector(t(CP)), each=m)
  e <- rnorm(nclust*Tp*m, mean=0, sd=sqrt(sig2e))
  
  # Combine terms to get outcome data for SW design
  Y <- betavec + Xvec*theta + CPvec + e
  Y0 <- betavec + Xvec*0 + CPvec + e
  
  # Create data frame with everything needed for fitting the model
  datSW <- data.frame(Y=Y, Y0=Y0, cluster=as.factor(clusterind), time=as.factor(perind), clustper=as.factor(clustperind), treat=Xvec)
  
  # Get subset of data corresponding to embedded basic staircase design
  datSC <- datSW[!is.na(XvecSC),]
  
  dats <- list(datSW, datSC)
  return(dats)
}

set.seed(3004)

# Calculate empirical power for the designs under exchangeable correlation
res20 <- emp_power(1000, SWdesmat(5,1), 8, 20, 0.03, 1, 0.15)
res30 <- emp_power(1000, SWdesmat(5,1), 8, 30, 0.03, 1, 0.15)

# Empirical relative efficiency
res20$emp_var_SW/res20$emp_var_SC
res20$emp_var_SW/res30$emp_var_SC

# Calculate empirical power for the designs under discrete-time decay correlation
res20_DTD <- emp_power(1000, SWdesmat(5,1), 8, 20, 0.032, 0.97, 0.15)
res30_DTD <- emp_power(1000, SWdesmat(5,1), 8, 30, 0.032, 0.97, 0.15)

# Empirical relative efficiency
res20_DTD$emp_var_SW/res20_DTD$emp_var_SC
res20_DTD$emp_var_SW/res30_DTD$emp_var_SC
