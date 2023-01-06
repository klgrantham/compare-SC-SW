# Comparison plots: staircase design vs. stepped wedge design
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

source('functions_releff.R')

library(tidyverse)
library(reshape2)

gridvals_small <- function(m_SW, S_SW, reps_SW, m_SC, S_SC, reps_SC,
                           pre_SC, post_SC, corrtype, pereff, fullrange=FALSE){
  # Compare variances of complete SW and staircase designs, for a range of
  # correlation parameters
  # Inputs:
  #  m_SW - number of subjects per cluster-period for SW design
  #  S_SW - number of unique treatment sequences for SW design
  #  reps_SW - number of times each sequence is repeated for SW design
  #  m_SC - number of subjects per cluster-period for SC design
  #  S_SC - number of unique treatment sequences for SC design
  #  reps_SC - number of times each sequence is repeated for SC design
  #  pre_SC - number of pre-switch measurement periods for SC design
  #  post_SC - number of post-switch measurement periods for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  # Output:
  #  Relative variances (vartheta_SC/vartheta_SW) and relative efficiencies
  #  (vartheta_SW/vartheta_SC) for a range of rho and r values
  
  if(fullrange==TRUE){
    rhovals <- seq(0.01, 0.99, 0.01)
    rvals <- seq(0.05, 1, 0.05)
  }else{
    rhovals <- seq(0.01, 0.25, 0.01)
    rvals <- seq(0.2, 1, 0.05)
  }
  vars <- expand.grid(rho=rhovals, r=rvals)
  
  vars$varSCmat <- with(
    vars,
    sapply(1:nrow(vars), function(j){
      CRTvartheta(m_SC, SCdesmat(S_SC, reps_SC, pre_SC, post_SC),
               rho[j], r[j], corrtype, pereff)
    }
    )
  )
  vars$varSW <- with(
    vars,
    sapply(1:nrow(vars), function(j){
      CRTvartheta(m_SW, SWdesmat(S_SW, reps_SW),
               rho[j], r[j], corrtype, pereff)
    }
    )
  )
  
  if(pereff=="cat"){
    vars <- vars %>%
      mutate(
        varSC = VarSCcat(m_SC, S_SC, reps_SC, rho, r),
      )
  }else if(pereff=="lin"){
    vars <- vars %>%
      mutate(
        varSC = VarSClin(m_SC, S_SC, reps_SC, rho, r),
      )
  }
  
  vars <- vars %>%
    mutate(
      relvarSCSW = varSC/varSW,
      relvarSCSWmat = varSCmat/varSW,
      releffSCSW = varSW/varSC
    )
  return(vars)
}

releffSCSW_grid_multiplot_corr <- function(S_SW, K_SW, S_SC, K_SC,
                                           pre_SC, post_SC, pereff,
                                           fixedscale=FALSE, limits=c(1,2),
                                           breaks=seq(1,2,0.5), genK=FALSE,
                                           fullrange=FALSE){
  # Compare variances of complete SW and staircase designs, for a range of
  # correlation parameters
  # Inputs:
  #  m_SW - number of subjects per cluster-period for SW design
  #  S_SW - number of unique treatment sequences for SW design
  #  K_SW - number of times each sequence is repeated for SW design
  #  m_SC - number of subjects per cluster-period for SC design
  #  S_SC - number of unique treatment sequences for SC design
  #  K_SC - number of times each sequence is repeated for SC design
  #  pre_SC - number of pre-switch measurement periods for SC design
  #  post_SC - number of post-switch measurement periods for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  # Output:
  #  Contour plot of relative variances (vartheta_SC/vartheta_SW)
  
  m1 <- 10
  m2 <- 100
  relvars_m1_BE <- gridvals_small(m1, S_SW, K_SW, m1, S_SC, K_SC,
                                  pre_SC, post_SC, 0, pereff, fullrange)
  relvars_m1_BE$m <- "m = 10"
  relvars_m1_BE$corrname <- "Block-exchangeable"
  relvars_m2_BE <- gridvals_small(m2, S_SW, K_SW, m2, S_SC, K_SC,
                                  pre_SC, post_SC, 0, pereff, fullrange)
  relvars_m2_BE$m <- "m = 100"
  relvars_m2_BE$corrname <- "Block-exchangeable"
  relvars_m1_DTD <- gridvals_small(m1, S_SW, K_SW, m1, S_SC, K_SC,
                                   pre_SC, post_SC, 1, pereff, fullrange)
  relvars_m1_DTD$m <- "m = 10"
  relvars_m1_DTD$corrname <- "Discrete-time decay"
  relvars_m2_DTD <- gridvals_small(m2, S_SW, K_SW, m2, S_SC, K_SC,
                                   pre_SC, post_SC, 1, pereff, fullrange)
  relvars_m2_DTD$m <- "m = 100"
  relvars_m2_DTD$corrname <- "Discrete-time decay"
  
  relvars <- bind_rows(
    relvars_m1_BE,
    relvars_m2_BE,
    relvars_m1_DTD,
    relvars_m2_DTD
  )
  
  if(fixedscale==TRUE){
    fillopt <- scale_fill_viridis_c(name="Relative efficiency ", direction=-1,
                                    limits=limits, breaks=breaks)
  }else{
    fillopt <- scale_fill_viridis_c(name="Relative efficiency ", direction=-1)
  }
  
  p <- ggplot(relvars, aes(x=r, y=rho)) +
    geom_tile(aes(fill=releffSCSW)) +
    fillopt +
    facet_grid(
      m ~ corrname
    ) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(aspect.ratio=3/8,
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=14), legend.text=element_text(size=14),
          legend.position="bottom",
          plot.title=element_text(hjust=0.5, size=16),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
          strip.background = element_rect(
            color="white", fill="white", linetype="solid"
          ),
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size=12)) +
    coord_fixed() + xlab(expression(paste("Cluster autocorrelation, ", r))) + ylab(expression(paste("Within-period ICC, ", rho))) +
    ggtitle(bquote(paste("Relative efficiency, ",
                         var(hat(theta))[paste("SW(", .(S_SW), ",", .(K_SW), ")")]/
                           var(hat(theta))[paste("SC(", .(S_SC), ",", .(K_SC), ",", .(pre_SC), ",", .(post_SC), ")")])))
  
  if(genK==TRUE){
    p <- p + ggtitle(bquote(paste("Relative efficiency, ",
                                  var(hat(theta))[paste("SW(", .(S_SW), ",k)")]/
                                    var(hat(theta))[paste("SC(", .(S_SC), ",k,", .(pre_SC), ",", .(post_SC), ")")])))
  }else{
    p <- p
  }
  
  rng <- ifelse(fullrange, "full", "restricted")

  ggsave(paste0("plots/releff_SC_", S_SC, K_SC, pre_SC, post_SC, "_vs_SW_",
                S_SW, K_SW, "_", pereff, "_", rng, ".jpg"),
         p, width=9, height=5, units="in", dpi=800)
  return(p)
}

releffSCSW_grid_multiplot_diffm <- function(S_SW, K_SW, S_SC, K_SC,
                                            pre_SC, post_SC, pereff,
                                            fixedscale=FALSE, limits=c(1,2),
                                            breaks=seq(1,2,0.5), genK=FALSE,
                                            fullrange=FALSE){
  # Compare variances of complete SW and staircase designs, for a range of
  # correlation parameters
  # Inputs:
  #  m_SW - number of subjects per cluster-period for SW design
  #  S_SW - number of unique treatment sequences for SW design
  #  K_SW - number of times each sequence is repeated for SW design
  #  m_SC - number of subjects per cluster-period for SC design
  #  S_SC - number of unique treatment sequences for SC design
  #  K_SC - number of times each sequence is repeated for SC design
  #  pre_SC - number of pre-switch measurement periods for SC design
  #  post_SC - number of post-switch measurement periods for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  # Output:
  #  Contour plot of relative variances (vartheta_SC/vartheta_SW)
  
  m1_SW <- 10
  m1_SC <- ((S_SW+1)/2)*m1_SW
  m2_SW <- 100
  m2_SC <- ((S_SW+1)/2)*m2_SW
  relvars_m1_BE <- gridvals_small(m1_SW, S_SW, K_SW, m1_SC, S_SC, K_SC,
                                  pre_SC, post_SC, 0, pereff, fullrange)
  relvars_m1_BE$m <- "m1"
  relvars_m1_BE$corrname <- "Block-exchangeable"
  relvars_m2_BE <- gridvals_small(m2_SW, S_SW, K_SW, m2_SC, S_SC, K_SC,
                                  pre_SC, post_SC, 0, pereff, fullrange)
  relvars_m2_BE$m <- "m2"
  relvars_m2_BE$corrname <- "Block-exchangeable"
  relvars_m1_DTD <- gridvals_small(m1_SW, S_SW, K_SW, m1_SC, S_SC, K_SC,
                                   pre_SC, post_SC, 1, pereff, fullrange)
  relvars_m1_DTD$m <- "m1"
  relvars_m1_DTD$corrname <- "Discrete-time decay"
  relvars_m2_DTD <- gridvals_small(m2_SW, S_SW, K_SW, m2_SC, S_SC, K_SC,
                                   pre_SC, post_SC, 1, pereff, fullrange)
  relvars_m2_DTD$m <- "m2"
  relvars_m2_DTD$corrname <- "Discrete-time decay"
  relvars <- bind_rows(
    relvars_m1_BE,
    relvars_m2_BE,
    relvars_m1_DTD,
    relvars_m2_DTD
  )
  
  if(fixedscale==TRUE){
    fillopt <- scale_fill_viridis_c(name="Relative efficiency ", direction=-1,
                                    limits=limits, breaks=breaks)
  }else{
    fillopt <- scale_fill_viridis_c(name="Relative efficiency ", direction=-1)
  }
  
  m.labs <- c(paste("mSC = ", m1_SC, "\nmSW = 10"), paste("mSC = ", m2_SC, "\nmSW = 100"))
  #  m.labs <- c("mSC=(S+1)mSW/2,\nmSW=10", "mSC=(S+1)mSW/2,\nmSW=100")
  names(m.labs) <- c("m1", "m2")
  
  p <- ggplot(relvars, aes(x=r, y=rho)) + 
    geom_tile(aes(fill=releffSCSW)) +
    facet_grid(
      m ~ corrname,
      labeller = labeller(m = m.labs)
    ) +
    fillopt +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(aspect.ratio=3/8,
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=14), legend.text=element_text(size=14),
          legend.position="bottom",
          plot.title=element_text(hjust=0.5, size=16),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
          strip.background = element_rect(
            color="white", fill="white", linetype="solid"
          ),
          strip.text.x = element_text(size = 10),
          strip.text.y = element_text(size=10)) +
    coord_fixed() + xlab(expression(paste("Cluster autocorrelation, ", r))) + ylab(expression(paste("Within-period ICC, ", rho))) +
    ggtitle(bquote(paste("Relative efficiency, ",
                         var(hat(theta))[paste("SW(", .(S_SW), ",", .(K_SW), ")")]/
                           var(hat(theta))[paste("SC(", .(S_SC), ",", .(K_SC), ",", .(pre_SC), ",", .(post_SC), ")")])))
  
    if(genK==TRUE){
      p <- p + ggtitle(bquote(paste("Relative efficiency, ",
                                    var(hat(theta))[paste("SW(", .(S_SW), ",k)")]/
                                      var(hat(theta))[paste("SC(", .(S_SC), ",k,", .(pre_SC), ",", .(post_SC), ")")])))
    }else{
      p <- p
    }
  
  rng <- ifelse(fullrange, "full", "restricted")
  
  ggsave(paste0("plots/releff_diffm_SC_", S_SC, K_SC, pre_SC, post_SC, "_vs_SW_",
                S_SW, K_SW, "_", pereff, "_", rng, ".jpg"),
         p, width=9, height=5, units="in", dpi=800)
  return(p)
}

## OLD ##
# gridvals_small <- function(m_SW, S_SW, reps_SW, m_SC, S_SC, reps_SC,
#                            pre_SC, post_SC, corrtype, pereff){
#   # Compare variances of complete SW and staircase designs, for a range of
#   # correlation parameters
#   # Inputs:
#   #  m_SW - number of subjects per cluster-period for SW design
#   #  S_SW - number of unique treatment sequences for SW design
#   #  reps_SW - number of times each sequence is repeated for SW design
#   #  m_SC - number of subjects per cluster-period for SC design
#   #  S_SC - number of unique treatment sequences for SC design
#   #  reps_SC - number of times each sequence is repeated for SC design
#   #  pre_SC - number of pre-switch measurement periods for SC design
#   #  post_SC - number of post-switch measurement periods for SC design
#   #  corrtype - within-cluster correlation structure type
#   #             (0=block-exchangeable, 1=exponential decay)
#   #  pereff - time period effect type
#   #           ('cat'=categorical period effects, 'lin'=linear period effects)
#   # Output:
#   #  Relative variances (vartheta_SC/vartheta_SW) for a range of rho and r values
#   
#   rho0seq <- seq(0.01, 0.25, 0.01)
#   rseq <- seq(0.2, 0.95, 0.05)
#   
#   SCSWvars <- matrix(data=NA, nrow=length(rho0seq), ncol=length(rseq))
#   for(i in 1:length(rho0seq)) {
#     for(rind in 1:length(rseq)) {
#       SCSWvars[i,rind] <-CRTvartheta(m_SC, SCdesmat(S_SC, reps_SC, pre_SC, post_SC),
#                                   rho0seq[i], rseq[rind], corrtype=corrtype, pereff=pereff)/
#         CRTvartheta(m_SW, SWdesmat(S_SW, reps_SW), rho0seq[i],
#                  rseq[rind], corrtype=corrtype, pereff=pereff)
#     }
#   }
#   
#   # Plot the results using a contour plot
#   SCSWvars<-round(SCSWvars, 2)
#   meltSCSWvars <- melt(SCSWvars)
#   
#   names(meltSCSWvars)[names(meltSCSWvars)=="Var1"] <- "rho"
#   names(meltSCSWvars)[names(meltSCSWvars)=="Var2"] <- "r"
#   
#   rhovec <- as.vector(matrix(data=rho0seq, nrow=length(rho0seq), ncol=length(rseq), byrow=FALSE))
#   rvec <- as.vector(matrix(data=rseq, nrow=length(rho0seq), ncol=length(rseq), byrow=TRUE))
#   meltSCSWvars$rhoseq <- rhovec
#   meltSCSWvars$rseq <- rvec
#   return(meltSCSWvars)
# }

## OLD ##
varSCSW_grid_small_plot <- function(m_SW, S_SW, reps_SW, m_SC, S_SC, reps_SC,
                                    pre_SC, post_SC, corrtype, pereff){
  # Compare variances of complete SW and staircase designs, for a range of
  # correlation parameters
  # Inputs:
  #  m_SW - number of subjects per cluster-period for SW design
  #  S_SW - number of unique treatment sequences for SW design
  #  reps_SW - number of times each sequence is repeated for SW design
  #  m_SC - number of subjects per cluster-period for SC design
  #  S_SC - number of unique treatment sequences for SC design
  #  reps_SC - number of times each sequence is repeated for SC design
  #  pre_SC - number of pre-switch measurement periods for SC design
  #  post_SC - number of post-switch measurement periods for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  # Output:
  #  Contour plot of relative variances (vartheta_SC/vartheta_SW)
  
  relvars <- gridvals_small(m_SW, S_SW, reps_SW, m_SC, S_SC, reps_SC, pre_SC, post_SC,
                            corrtype, pereff)  
  
  p <- ggplot(relvars, aes(x=rseq, y=rhoseq)) + 
    geom_tile(aes(fill=value)) + 
    scale_fill_gradientn(colours=c("yellow","red")) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(aspect.ratio=3/8, legend.position="bottom", legend.key.size=unit(1, "cm"), 
          legend.text=element_text(size=12), 
          legend.background = element_rect(fill="grey95"),
          plot.title=element_text(hjust=0.5, size=16),
          axis.title=element_text(size=14), axis.text=element_text(size=14)) +
    coord_fixed() + xlab(expression(paste("Cluster autocorrelation, ", r))) + ylab(expression(paste("Within-period ICC, ", rho))) +
    geom_text(aes(rseq, rhoseq, label=round(value,2)), color="black", size=3)  +
    ggtitle(bquote(paste("Relative variance of treatment effect estimators, ",
                         Var(hat(theta))[paste("SC(", .(S_SC), ",", .(reps_SC), ",", .(pre_SC), ",", .(post_SC), ")")]/
                           Var(hat(theta))[paste("SW(", .(S_SW), ",", .(reps_SW), ")")])))
  corrname <- ifelse(corrtype==0, "BE", "DTD")
  ggsave(paste0("plots/SC_", S_SC, reps_SC, pre_SC, post_SC, "_vs_SW_",
                S_SW, reps_SW, "_", corrname, "_", pereff, ".jpg"),
         p, width=9, height=7, units="in", dpi=800)
  return(p)
}

## OLD ##
varSCSW_grid_multiplot <- function(S_SW, reps_SW, S_SC, reps_SC,
                                   pre_SC, post_SC, corrtype){
  # Compare variances of complete SW and staircase designs, for a range of
  # correlation parameters
  # Inputs:
  #  m_SW - number of subjects per cluster-period for SW design
  #  S_SW - number of unique treatment sequences for SW design
  #  reps_SW - number of times each sequence is repeated for SW design
  #  m_SC - number of subjects per cluster-period for SC design
  #  S_SC - number of unique treatment sequences for SC design
  #  reps_SC - number of times each sequence is repeated for SC design
  #  pre_SC - number of pre-switch measurement periods for SC design
  #  post_SC - number of post-switch measurement periods for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  # Output:
  #  Contour plot of relative variances (vartheta_SC/vartheta_SW)
  
  m1 <- 10
  m2 <- 100
  relvars_m1_cat <- gridvals_small(m1, S_SW, reps_SW, m1, S_SC, reps_SC,
                                   pre_SC, post_SC, corrtype, 'cat')
  relvars_m1_cat$m <- "m = 10"
  relvars_m1_cat$pereff <- "Categorical time"
  relvars_m2_cat <- gridvals_small(m2, S_SW, reps_SW, m2, S_SC, reps_SC,
                                   pre_SC, post_SC, corrtype, 'cat')
  relvars_m2_cat$m <- "m = 100"
  relvars_m2_cat$pereff <- "Categorical time"
  relvars_m1_lin <- gridvals_small(m1, S_SW, reps_SW, m1, S_SC, reps_SC,
                                   pre_SC, post_SC, corrtype, 'lin')
  relvars_m1_lin$m <- "m = 10"
  relvars_m1_lin$pereff <- "Linear time"
  relvars_m2_lin <- gridvals_small(m2, S_SW, reps_SW, m2, S_SC, reps_SC,
                                   pre_SC, post_SC, corrtype, 'lin')
  relvars_m2_lin$m <- "m = 100"
  relvars_m2_lin$pereff <- "Linear time"
  relvars <- bind_rows(
    relvars_m1_cat,
    relvars_m2_cat,
    relvars_m1_lin,
    relvars_m2_lin
  )
  
  p <- ggplot(relvars, aes(x=rseq, y=rhoseq)) + 
    geom_tile(aes(fill=value)) +
    facet_grid(
      m ~ pereff
    ) +
    scale_fill_viridis_c(name="Relative variance", direction=-1) +
    #    scale_fill_gradientn(colours=c("yellow","red")) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(aspect.ratio=3/8,
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=14), legend.text=element_text(size=14),
          legend.position="bottom",
          plot.title=element_text(hjust=0.5, size=16),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
          strip.background = element_rect(
            color="white", fill="white", linetype="solid"
          ),
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size=12)) +
    coord_fixed() + xlab(expression(paste("Cluster autocorrelation, ", r))) + ylab(expression(paste("Within-period ICC, ", rho))) +
    #    geom_text(aes(rseq, rhoseq, label=round(value,2)), color="black", size=3)  +
    ggtitle(bquote(paste("Relative variance of treatment effect estimators, ",
                         Var(hat(theta))[paste("SC(", .(S_SC), ",", .(reps_SC), ",", .(pre_SC), ",", .(post_SC), ")")]/
                           Var(hat(theta))[paste("SW(", .(S_SW), ",", .(reps_SW), ")")])))
  corrname <- ifelse(corrtype==0, "BE", "DTD")
  ggsave(paste0("plots/multiplot_SC_", S_SC, reps_SC, pre_SC, post_SC, "_vs_SW_",
                S_SW, reps_SW, "_", corrname, ".jpg"),
         p, width=9, height=5, units="in", dpi=800)
  return(p)
}

## OLD ##
gridvals_varSC <- function(m_SC, S_SC, reps_SC, pre_SC, post_SC, corrtype, pereff){
  # Get variances of staircase design, for a range of correlation parameters
  # Inputs:
  #  m_SC - number of subjects per cluster-period for SC design
  #  S_SC - number of unique treatment sequences for SC design
  #  reps_SC - number of times each sequence is repeated for SC design
  #  pre_SC - number of pre-switch measurement periods for SC design
  #  post_SC - number of post-switch measurement periods for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  # Output:
  #  Variances for a range of rho and r values
  
  rho0seq <- seq(0.05, 0.95, 0.05)
  rseq <- seq(0.05, 0.95, 0.05)
  
  SCvars <- matrix(data=NA, nrow=length(rho0seq), ncol=length(rseq))
  for(i in 1:length(rho0seq)) {
    for(rind in 1:length(rseq)) {
      SCvars[i,rind] <-CRTvartheta(m_SC, SCdesmat(S_SC, reps_SC, pre_SC, post_SC),
                                rho0seq[i], rseq[rind], corrtype=corrtype, pereff=pereff)
    }
  }
  
  # Plot the results using a contour plot
  SCvars<-round(SCvars, 4)
  meltSCvars <- melt(SCvars)
  
  names(meltSCvars)[names(meltSCvars)=="Var1"] <- "rho"
  names(meltSCvars)[names(meltSCvars)=="Var2"] <- "r"
  
  rhovec <- as.vector(matrix(data=rho0seq, nrow=length(rho0seq), ncol=length(rseq), byrow=FALSE))
  rvec <- as.vector(matrix(data=rseq, nrow=length(rho0seq), ncol=length(rseq), byrow=TRUE))
  meltSCvars$rhoseq <- rhovec
  meltSCvars$rseq <- rvec
  return(meltSCvars)
}

## OLD ##
varSC_grid_multiplot <- function(S1_SC, S2_SC, reps_SC, pre_SC, post_SC, corrtype, pereff){
  # Display variances of staircase designs, for a range of
  # correlation parameters
  # Inputs:
  #  S1_SC - number of unique treatment sequences for SC design (left column)
  #  S2_SC - number of unique treatment sequences for SC design (right column)
  #  reps_SC - number of times each sequence is repeated for SC design
  #  pre_SC - number of pre-switch measurement periods for SC design
  #  post_SC - number of post-switch measurement periods for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  # Output:
  #  Contour plot of relative variances (vartheta_SC/vartheta_SW)
  
  m1 <- 10
  m2 <- 100
  S1 <- paste("S = ", S1_SC)
  S2 <- paste("S = ", S2_SC)
  
  vars_m1_S1 <- gridvals_varSC(m1, S1_SC, reps_SC, pre_SC, post_SC, corrtype, pereff)
  vars_m1_S1$m <- "m = 10"
  vars_m1_S1$Sname <- S1
  vars_m2_S1 <- gridvals_varSC(m2, S1_SC, reps_SC, pre_SC, post_SC, corrtype, pereff)
  vars_m2_S1$m <- "m = 100"
  vars_m2_S1$Sname <- S1
  vars_m1_S2 <- gridvals_varSC(m1, S2_SC, reps_SC, pre_SC, post_SC, corrtype, pereff)
  vars_m1_S2$m <- "m = 10"
  vars_m1_S2$Sname <- S2
  vars_m2_S2 <- gridvals_varSC(m2, S2_SC, reps_SC, pre_SC, post_SC, corrtype, pereff)
  vars_m2_S2$m <- "m = 100"
  vars_m2_S2$Sname <- S2
  vars <- bind_rows(
    vars_m1_S1,
    vars_m2_S1,
    vars_m1_S2,
    vars_m2_S2
  )
  
  p <- ggplot(vars, aes(x=rseq, y=rhoseq)) + 
    geom_tile(aes(fill=value)) +
    facet_grid(
      m ~ Sname
    ) +
    scale_fill_viridis_c(name="Variance", direction=-1) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(aspect.ratio=3/8,
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=14), legend.text=element_text(size=14),
          legend.position="bottom",
          plot.title=element_text(hjust=0.5, size=16),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
          strip.background = element_rect(
            color="white", fill="white", linetype="solid"
          ),
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size=12)) +
    coord_fixed() + xlab(expression(paste("Cluster autocorrelation, ", r))) + ylab(expression(paste("Within-period ICC, ", rho))) +
    ggtitle(bquote(paste("Variance of treatment effect estimator, ",
                         Var(hat(theta))[paste("SC(S,", .(reps_SC), ",", .(pre_SC), ",", .(post_SC), ")")])))
  ggsave(paste0("plots/multiplot_SC_", S1_SC, "vs", S2_SC, reps_SC, pre_SC, post_SC, "_", corrtype, "_", pereff, ".jpg"),
         p, width=9, height=5, units="in", dpi=800)
  return(p)
}

## OLD ##
gridvals <- function(m_SW, S_SW, reps_SW, m_SC, S_SC, reps_SC,
                     pre_SC, post_SC, corrtype, pereff){
  # Compare variances of complete SW and staircase designs, for a range of
  # correlation parameters
  # Inputs:
  #  m_SW - number of subjects per cluster-period for SW design
  #  S_SW - number of unique treatment sequences for SW design
  #  reps_SW - number of times each sequence is repeated for SW design
  #  m_SC - number of subjects per cluster-period for SC design
  #  S_SC - number of unique treatment sequences for SC design
  #  reps_SC - number of times each sequence is repeated for SC design
  #  pre_SC - number of pre-switch measurement periods for SC design
  #  post_SC - number of post-switch measurement periods for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  # Output:
  #  Relative variances (vartheta_SC/vartheta_SW) for a range of rho and r values
  
  rho0seq <- seq(0.01, 0.99, 0.05)
  rseq <- seq(0.0, 0.95, 0.05)
  
  SCSWvars <- matrix(data=NA, nrow=length(rho0seq), ncol=length(rseq))
  for(i in 1:length(rho0seq)) {
    for(rind in 1:length(rseq)) {
      SCSWvars[i,rind] <-CRTvartheta(m_SC, SCdesmat(S_SC, reps_SC, pre_SC, post_SC),
                                  rho0seq[i], rseq[rind], corrtype=corrtype, pereff=pereff)/
        CRTvartheta(m_SW, SWdesmat(S_SW, reps_SW), rho0seq[i],
                 rseq[rind], corrtype=corrtype, pereff=pereff)
    }
  }
  
  # Plot the results using a contour plot
  SCSWvars<-round(SCSWvars, 2)
  meltSCSWvars <- melt(SCSWvars)
  
  names(meltSCSWvars)[names(meltSCSWvars)=="Var1"] <- "rho"
  names(meltSCSWvars)[names(meltSCSWvars)=="Var2"] <- "r"
  
  rhovec <- as.vector(matrix(data=rho0seq, nrow=length(rho0seq), ncol=length(rseq), byrow=FALSE))
  rvec <- as.vector(matrix(data=rseq, nrow=length(rho0seq), ncol=length(rseq), byrow=TRUE))
  meltSCSWvars$rhoseq <- rhovec
  meltSCSWvars$rseq <- rvec
  return(meltSCSWvars)
}

## OLD ##
varSCSW_grid_plot <- function(m_SW, S_SW, reps_SW, m_SC, S_SC, reps_SC,
                              pre_SC, post_SC, corrtype, pereff){
  # Compare variances of complete SW and staircase designs, for a range of
  # correlation parameters
  # Inputs:
  #  m_SW - number of subjects per cluster-period for SW design
  #  S_SW - number of unique treatment sequences for SW design
  #  reps_SW - number of times each sequence is repeated for SW design
  #  m_SC - number of subjects per cluster-period for SC design
  #  S_SC - number of unique treatment sequences for SC design
  #  reps_SC - number of times each sequence is repeated for SC design
  #  pre_SC - number of pre-switch measurement periods for SC design
  #  post_SC - number of post-switch measurement periods for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  # Output:
  #  Contour plot of relative variances (vartheta_SC/vartheta_SW)
  
  relvars <- gridvals(m_SW, S_SW, reps_SW, m_SC, S_SC, reps_SC, pre_SC, post_SC,
                      corrtype, pereff)  
  
  myplot <- ggplot(relvars, aes(x=rseq, y=rhoseq)) + 
    geom_tile(aes(fill=value)) + 
    scale_fill_gradientn(colours=c("yellow","red")) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(aspect.ratio=3/8, legend.position="none", legend.key.size=unit(1, "cm"), 
          legend.text=element_text(size=12), 
          legend.background = element_rect(fill="grey95")) +
    coord_fixed() + xlab("Cluster autocorrelation, r") +  ylab("Within-period ICC") +
    geom_text(aes(rseq, rhoseq, label=round(value,2)), color="black", size=3)
  
  return(myplot)
}

## OLD ##
varSCSW_line_plot <- function(m_SW, S_SW, reps_SW, m_SC, S_SC, reps_SC,
                              pre_SC, post_SC, corrtype, pereff, title=""){
  
  rhovals <- seq(0.01, 0.2, 0.005)
  rvals <- c(0.25, 0.5, 0.75, 0.95, 1.0)
  relvars <- expand.grid(rho=rhovals, r=rvals)
  relvars$varSC <- with(
    relvars,
    sapply(1:nrow(relvars), function(j){
      CRTvartheta(m_SC, SCdesmat(S_SC, reps_SC, pre_SC, post_SC),
               rho[j], r[j], corrtype, pereff)
    }
    )
  )
  relvars$varSW <- with(
    relvars,
    sapply(1:nrow(relvars), function(j){
      CRTvartheta(m_SW, SWdesmat(S_SW, reps_SW),
               rho[j], r[j], corrtype, pereff)
    }
    )
  )
  relvars <- relvars %>%
    mutate(relvarSCSW = varSC/varSW,
           rfac=as.factor(r)) %>%
    select(c(rho, rfac, relvarSCSW))
  
  p <- ggplot(data=relvars, aes(x=rho, y=relvarSCSW, colour=rfac)) +
    geom_line(size=1.2) +
    geom_hline(yintercept=1, linetype="dashed") +
    #      expand_limits(y=ylimits) +
    xlab("Within-period ICC") +
    ylab("Relative variance") +
    labs(title=title, colour="Cluster autocorrelation") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=12),
          axis.title=element_text(size=10), axis.text=element_text(size=10),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=12), legend.text=element_text(size=12),
          legend.position="bottom")
  return(p)
}

## OLD ##
varSCSW_multi_plot <- function(S_SW, reps_SW, S_SC, reps_SC, pre_SC, post_SC, corrtype){
  # Compare variances of complete SW and staircase designs, for a range of
  # correlation parameters
  # Inputs:
  #  S_SW - number of unique treatment sequences for SW design
  #  reps_SW - number of times each sequence is repeated for SW design
  #  S_SC - number of unique treatment sequences for SC design
  #  reps_SC - number of times each sequence is repeated for SC design
  #  pre_SC - number of pre-switch measurement periods for SC design
  #  post_SC - number of post-switch measurement periods for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  # Output:
  #  Multiplot of relative variances (vartheta_SC/vartheta_SW), for m=10 and 100
  #  and categorical and linear period effects
  
  m1 <- 10
  m2 <- 100
  
  # m=10, categorical period effects
  p1 <- varSCSW_line_plot(
    m1, S_SW, reps_SW,
    m1, S_SC, reps_SC, pre_SC, post_SC,
    corrtype, 'cat',
    bquote(paste("m = ", .(m1), ", ", "categorical period effects"))
  )
  # m=100, categorical period effects
  p2 <- varSCSW_line_plot(
    m2, S_SW, reps_SW,
    m2, S_SC, reps_SC, pre_SC, post_SC,
    corrtype, 'cat',
    bquote(paste("m = ", .(m2), ", ", "categorical period effects"))
  )
  # m=10, linear period effects
  p3 <- varSCSW_line_plot(
    m1, S_SW, reps_SW,
    m1, S_SC, reps_SC, pre_SC, post_SC,
    corrtype, 'lin',
    bquote(paste("m = ", .(m1), ", ", "linear period effects"))
  )
  # m=100, linear period effects
  p4 <- varSCSW_line_plot(
    m2, S_SW, reps_SW,
    m2, S_SC, reps_SC, pre_SC, post_SC,
    corrtype, 'lin',
    bquote(paste("m = ", .(m2), ", ", "linear period effects"))
  )
  mylegend <- g_legend(p1)
  title <- bquote(paste("Relative variance of treatment effect estimators, ",
                        Var(hat(theta))[paste("SC(", .(S_SC), ",", .(reps_SC), ",", .(pre_SC), ",", .(post_SC), ")")]/
                          Var(hat(theta))[paste("SW(", .(S_SW), ",", .(reps_SW), ")")]))
  p1to4 <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title=title)
  corrname <- ifelse(corrtype==0, "BE", "DTD")
  ggsave(paste0("plots/SC", S_SC, reps_SC, pre_SC, post_SC, "_vs_SW",
                S_SW, reps_SW, "_", corrname, ".jpg"),
         p1to4, width=9, height=7, units="in", dpi=800)
  
  return(p1to4)
}

## OLD ##
varSCSW_line_plot_sequences <- function(S_vals, r, m_SW, reps_SW, m_SC, reps_SC,
                                        pre_SC, post_SC, corrtype, pereff, title=""){
  
  rhovals <- seq(0.01, 0.2, 0.005)
  r <- r
  relvars <- expand.grid(S=S_vals, rho=rhovals, r=r)
  relvars$varSC <- with(
    relvars,
    sapply(1:nrow(relvars), function(j){
      CRTvartheta(m_SC, SCdesmat(S[j], reps_SC, pre_SC, post_SC),
               rho[j], r[j], corrtype, pereff)
    }
    )
  )
  relvars$varSW <- with(
    relvars,
    sapply(1:nrow(relvars), function(j){
      CRTvartheta(m_SW, SWdesmat(S[j], reps_SW),
               rho[j], r[j], corrtype, pereff)
    }
    )
  )
  relvars <- relvars %>%
    mutate(relvarSCSW = varSC/varSW,
           Sfac=as.factor(S)) %>%
    select(c(rho, Sfac, relvarSCSW))
  
  p <- ggplot(data=relvars, aes(x=rho, y=relvarSCSW, colour=Sfac)) +
    geom_line(size=1.2) +
    geom_hline(yintercept=1, linetype="dashed") +
    #      expand_limits(y=ylimits) +
    xlab("Within-period ICC") +
    ylab("Relative variance") +
    labs(title=title, colour="Sequences") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=12),
          axis.title=element_text(size=10), axis.text=element_text(size=10),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=12), legend.text=element_text(size=12),
          legend.position="bottom")
  return(p)
}

## OLD ##
varSCSW_multi_plot_sequences <- function(S_vals, r, reps_SW, reps_SC, pre_SC, post_SC, corrtype){
  # Compare variances of complete SW and staircase designs, for a range of
  # correlation parameters
  # Inputs:
  #  S_vals - numbers of unique treatment sequences, to form lines on plots, e.g. c(3, 5, 10, 20)
  #  r - cluster autocorrelation value for all scenarios
  #  reps_SW - number of times each sequence is repeated for SW design
  #  reps_SC - number of times each sequence is repeated for SC design
  #  pre_SC - number of pre-switch measurement periods for SC design
  #  post_SC - number of post-switch measurement periods for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  # Output:
  #  Multiplot of relative variances (vartheta_SC/vartheta_SW), for m=10 and 100
  #  and categorical and linear period effects
  
  m1 <- 10
  m2 <- 100
  
  # m=10, categorical period effects
  p1 <- varSCSW_line_plot_sequences(
    S_vals, r,
    m1, reps_SW,
    m1, reps_SC, pre_SC, post_SC,
    corrtype, 'cat',
    bquote(paste("m = ", .(m1), ", ", "categorical period effects"))
  )
  # m=100, categorical period effects
  p2 <- varSCSW_line_plot_sequences(
    S_vals, r,
    m2, reps_SW,
    m2, reps_SC, pre_SC, post_SC,
    corrtype, 'cat',
    bquote(paste("m = ", .(m2), ", ", "categorical period effects"))
  )
  # m=10, linear period effects
  p3 <- varSCSW_line_plot_sequences(
    S_vals, r, 
    m1, reps_SW,
    m1, reps_SC, pre_SC, post_SC,
    corrtype, 'lin',
    bquote(paste("m = ", .(m1), ", ", "linear period effects"))
  )
  # m=100, linear period effects
  p4 <- varSCSW_line_plot_sequences(
    S_vals, r,
    m2, reps_SW,
    m2, reps_SC, pre_SC, post_SC,
    corrtype, 'lin',
    bquote(paste("m = ", .(m2), ", ", "linear period effects"))
  )
  mylegend <- g_legend(p1)
  title <- bquote(paste("Relative variance of treatment effect estimators, ",
                        Var(hat(theta))[paste("SC(S,", .(reps_SC), ",", .(pre_SC), ",", .(post_SC), ")")]/
                          Var(hat(theta))[paste("SW(S,", .(reps_SW), ")")]))
  p1to4 <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title=title)
  corrname <- ifelse(corrtype==0, "BE", "DTD")
  ggsave(paste0("plots/SC_S", reps_SC, pre_SC, post_SC, "_vs_SW_S",
                reps_SW, "_", corrname, ".jpg"),
         p1to4, width=9, height=7, units="in", dpi=800)
  return(p1to4)
}

VarSW_alt <- function(m, S, K, rho0, r){
  totalvar <- 1
  sig2CP <- rho0*totalvar
  sig2E <- totalvar - sig2CP
  sig2 <- sig2E/m
  a <- sig2 + sig2CP
  b <- r*sig2CP
  (12*(a - b)*(a + S*b))/(K*(S^2 - 1)*(2*a + S*b))
}


## Relative efficiency, limited parameter range

# Embedded staircase vs stepped wedge
releffSCSW_grid_multiplot_corr(3, 1, 3, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2), genK=TRUE)
releffSCSW_grid_multiplot_corr(3, 1, 3, 1, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2), genK=TRUE)

releffSCSW_grid_multiplot_corr(10, 1, 10, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2), genK=TRUE)
releffSCSW_grid_multiplot_corr(10, 1, 10, 1, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2), genK=TRUE)

# Staircase with larger cluster-period size vs stepped wedge
releffSCSW_grid_multiplot_diffm(3, 1, 3, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5), genK=TRUE)
releffSCSW_grid_multiplot_diffm(3, 1, 3, 1, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5), genK=TRUE)

releffSCSW_grid_multiplot_diffm(10, 1, 10, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5), genK=TRUE)
releffSCSW_grid_multiplot_diffm(10, 1, 10, 1, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5), genK=TRUE)

# Extended staircase vs stepped wedge
releffSCSW_grid_multiplot_corr(3, 1, 3, 2, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.5))
releffSCSW_grid_multiplot_corr(3, 1, 3, 2, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.5))


## Relative efficiency, full parameter range

# Embedded staircase vs stepped wedge
releffSCSW_grid_multiplot_corr(3, 1, 3, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2), genK=TRUE, fullrange=TRUE)
releffSCSW_grid_multiplot_corr(3, 1, 3, 1, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2), genK=TRUE, fullrange=TRUE)

releffSCSW_grid_multiplot_corr(10, 1, 10, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2), genK=TRUE, fullrange=TRUE)
releffSCSW_grid_multiplot_corr(10, 1, 10, 1, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2), genK=TRUE, fullrange=TRUE)

# Staircase with larger cluster-period size vs stepped wedge
releffSCSW_grid_multiplot_diffm(3, 1, 3, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5), genK=TRUE, fullrange=TRUE)
releffSCSW_grid_multiplot_diffm(3, 1, 3, 1, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5), genK=TRUE, fullrange=TRUE)

releffSCSW_grid_multiplot_diffm(10, 1, 10, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5), genK=TRUE, fullrange=TRUE)
releffSCSW_grid_multiplot_diffm(10, 1, 10, 1, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5), genK=TRUE, fullrange=TRUE)

# Extended staircase vs stepped wedge
releffSCSW_grid_multiplot_corr(3, 1, 3, 2, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.5), fullrange=TRUE)
releffSCSW_grid_multiplot_corr(3, 1, 3, 2, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.5), fullrange=TRUE)


## Trial examples

# PROMPT trial
SWvar <- CRTvartheta(20, SWdesmat(5, 1), 0.03, 1, 0, 'cat')
SCvar_m20 <- VarSCcat(20, 5, 1, 0.03, 1)
releff_SW_SCm20 <- SWvar/SCvar_m20
pow(SCvar_m20, 0.4)

SCvar_m33 <- VarSCcat(33, 5, 1, 0.03, 1)
releff_SW_SCm33 <- SWvar/SCvar_m33
pow(SCvar_m33, 0.4)

SWvar_dtd <- CRTvartheta(20, SWdesmat(5, 1), 0.033, 0.951, 1, 'cat')
SCvar_m20_dtd <- VarSCcat(20, 5, 1, 0.033, 0.951)
releff_SW_SCm20_dtd <- SWvar_dtd/SCvar_m20_dtd
pow(SCvar_m20_dtd, 0.4)

SCvar_m33_dtd <- VarSCcat(33, 5, 1, 0.033, 0.951)
releff_SW_SCm33_dtd <- SWvar_dtd/SCvar_m33_dtd
pow(SCvar_m33_dtd, 0.4)

# INSPIRED trial
SWvar <- CRTvartheta(5, SWdesmat(6, 2), 0.05, 1, 0, 'cat')
pow(SWvar, 0.5)
SCvar_m5 <- VarSCcat(5, 6, 2, 0.05, 1)
releff_SW_SCm5 <- SWvar/SCvar_m5
pow(SCvar_m5, 0.5)

SCvar_m9 <- VarSCcat(9, 6, 2, 0.05, 1)
releff_SW_SCm9 <- SWvar/SCvar_m9
pow(SCvar_m9, 0.5)

SWvar_dtd <- CRTvartheta(5, SWdesmat(6, 2), 0.056, 0.952, 1, 'cat')
SCvar_m5_dtd <- VarSCcat(5, 6, 2, 0.056, 0.952)
releff_SW_SCm5_dtd <- SWvar_dtd/SCvar_m5_dtd
pow(SCvar_m5_dtd, 0.5)

SCvar_m9_dtd <- VarSCcat(9, 6, 2, 0.056, 0.952)
releff_SW_SCm9_dtd <- SWvar_dtd/SCvar_m9_dtd
pow(SCvar_m9_dtd, 0.5)
