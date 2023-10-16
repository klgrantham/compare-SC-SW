# Comparison plots: staircase design vs. stepped wedge design
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

source('functions_releff.R')

library(tidyverse)
library(reshape2)
library(viridis)
library(ggplot2)
library(scales)
library(colorspace)

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
#    fillopt <- scale_fill_viridis_c(name="Relative efficiency ", direction=-1,
#                                    limits=limits, breaks=breaks)
#    fillopt <- scale_fill_gradient2(name="Relative efficiency ", midpoint=1,
#                                    limits=limits, breaks=breaks)
#    fillopt <- scale_fill_gradientn(name="Relative efficiency ",
#                                    colours=c(muted("red"), "yellow",
#                                              "white",
#                                              "green", muted("blue")),
#                                    limits=limits, breaks=breaks)
#    fillopt <- scale_fill_distiller(name="Relative efficiency ",
#                                    type="div", palette="RdBu", direction=1, #palette="PuOr"
#                                    limits=limits, breaks=breaks)
    fillopt <- scale_fill_continuous_divergingx(name="Relative efficiency ",
                                                palette='RdBu', mid=1, rev=FALSE,
                                                limits=limits, breaks=breaks) #palette='RdBu' 'Temps' 'RdYlBu'
  }else{
#    fillopt <- scale_fill_viridis_c(name="Relative efficiency ", direction=-1)
#    fillopt <- scale_fill_gradient2(name="Relative efficiency ", midpoint=1)
#    fillopt <- scale_fill_gradientn(name="Relative efficiency ",
#                                    colours=c(muted("red"), "yellow",
#                                              "white",
#                                              "green", muted("blue")))
#    fillopt <- scale_fill_distiller(name="Relative efficiency ",
#                                    type="div", palette="RdBu", direction=1)
    fillopt <- scale_fill_continuous_divergingx(name="Relative efficiency ",
                                                palette='RdBu', mid=1, rev=FALSE)
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
#    fillopt <- scale_fill_viridis_c(name="Relative efficiency ", direction=-1,
#                                    limits=limits, breaks=breaks)
    fillopt <- scale_fill_continuous_divergingx(name="Relative efficiency ",
                                                palette='RdBu', mid=1, rev=FALSE,
                                                limits=limits, breaks=breaks)
  }else{
#    fillopt <- scale_fill_viridis_c(name="Relative efficiency ", direction=-1)
    fillopt <- scale_fill_continuous_divergingx(name="Relative efficiency ",
                                                palette='RdBu', mid=1, rev=FALSE)
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

corr_grid_multiplot <- function(fixedscale=TRUE, limits=c(0,1),
                                breaks=seq(0,1,0.2), fullrange=FALSE){
  
  if(fullrange==TRUE){
    rhovals <- seq(0.01, 0.99, 0.01)
    rvals <- seq(0.05, 1, 0.05)
    mvals <- c(10, 100)
  }else{
    rhovals <- seq(0.01, 0.25, 0.01)
    rvals <- seq(0.2, 1, 0.05)
    mvals <- c(10, 100)
  }
  vars <- expand.grid(rho=rhovals, r=rvals, mvals=mvals)
  
  vars$m <- ifelse(vars$mvals==10, "m = 10", "m = 100")
  
  vars$corr <- with(
    vars,
    sapply(1:nrow(vars), function(j){
      psi_corr(mvals[j], rho[j], r[j])
    }
    )
  )
  
  if(fixedscale==TRUE){
#    fillopt <- scale_fill_viridis(name="Correlation ", direction=-1,
#                                    limits=limits, breaks=breaks)
    fillopt <- scale_fill_distiller(name="Correlation", 
                                    palette="YlOrRd", direction=1,
                                    limits=limits, breaks=breaks)
  }else{
#    fillopt <- scale_fill_viridis(name="Correlation ", direction=-1)
    fillopt <- scale_fill_distiller(name="Correlation",
                                    palette="YlOrRd", direction=1)
  }

  p <- ggplot(vars, aes(x=r, y=rho)) + 
    geom_tile(aes(fill=corr)) +
    facet_grid(m ~ .) +
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
#    ) +
    coord_fixed() + xlab(expression(paste("Cluster autocorrelation, ", r))) + ylab(expression(paste("Within-period ICC, ", rho))) +
    ggtitle(expression(paste("Correlation between cluster-period means, ", psi)))
  
  rng <- ifelse(fullrange, "full", "restricted")
  
  ggsave(paste0("plots/corr_m10_m100_", rng, ".jpg"),
         p, width=6, height=5, units="in", dpi=800)
  return(p)
}

weight_grid_multiplot <- function(){
  
  psivals <- seq(0, 1, 0.01)
  Svals <- seq(3, 20, 1) # T = S+1, T-1 = S

  vals <- expand.grid(psi=psivals, S=Svals)
  
  vals$b <- with(
    vals,
    sapply(1:nrow(vals), function(j){
      ((S[j]+1)*psi[j])/(1 + (psi[j]*S[j]))
    }
    )
  )
  
  fillopt <- scale_fill_gradient2(name="Weight ", midpoint=0.5,
                                    limits=c(0,1), breaks=seq(0,1,0.25))
  fillopt <- scale_fill_continuous_divergingx(name="Weight ",
                                              palette='PRGn', mid=0.5, rev=TRUE,
                                              limits=c(0,1), breaks=seq(0,1,0.25)) #palette='RdBu' 'Temps' 'RdYlBu'

  p <- ggplot(vals, aes(x=psi, y=S)) + 
    geom_tile(aes(fill=b)) +
    fillopt +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(aspect.ratio=3/8,
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=14), legend.text=element_text(size=14),
          legend.position="bottom",
          plot.title=element_text(hjust=0.5, size=16),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
    ) +
    coord_fixed() + xlab(expression(paste("Correlation between cluster-period means, ", psi))) +
    ylab("Number of sequences, S") +
    ggtitle("Weight on horizontal comparisons in a stepped wedge design")
  
  ggsave(paste0("plots/weight_horizontal.jpg"),
         p, width=7, height=4, units="in", dpi=800)
  return(p)
}

VarSW_no_a <- function(S, K, psi){
  (12*(1 - psi)*(1 + S*psi))/(K*(S^2 - 1)*(2 + S*psi))
}

VarSCcat_no_a <- function(S, K, psi){
  #  a <- (1 + (m-1)*rho0)/m
  #  b <- r*rho0
  
  fracterm <- ((1 + sqrt(1-psi^2))^S - psi^S)/((1 + sqrt(1-psi^2))^S + psi^S)
  vartheta <- (2*(1-psi)^2)/(K*(S*(1-psi)-sqrt(1-psi^2)*fracterm))
  return(vartheta)
}

VarSClin_no_a <- function(S, K, psi){
  vartheta <- 2*((S^2+2) - (S^2-4)*psi)/(K*S*(S^2-1))
  return(vartheta)
}

Varvals_psi <- function(S, K){

  psivals <- seq(0, 1, 0.01)
  vars <- data.frame(psi=psivals)

  vars <- vars %>%
    mutate(
      SCcat = VarSCcat_no_a(S, K, psi),
      SClin = VarSClin_no_a(S, K, psi),
      SW = VarSW_no_a(S, K, psi),
      S = S
    )
  
  vars_long <- vars %>%
    pivot_longer(
      cols = c("SCcat", "SClin", "SW"),
      names_to = "design",
      values_to = "variance"
    )
  
  return(vars_long)  
}

Varvals_a_psi <- function(m_SW, m_SC, rho0, S, K){
  
  rvals <- seq(0, 1, 0.01)
  vars <- data.frame(r=rvals)
  
  vars <- vars %>%
    mutate(
      psiSW = psi_corr(m_SW, rho0, r),
      psiSC = psi_corr(m_SC, rho0, r),
      SCcat = VarSCcat(m_SC, S, K, rho0, r),
      SClin = VarSClin(m_SC, S, K, rho0, r),
      SW = VarSW_alt(m_SW, S, K, rho0, r),
      mSW = m_SW,
      mSC = m_SC,
      mlab = paste("mSW = ", mSW, "\nmSC = ", mSC),
      rho0 = rho0,
      S = S
    )
  
  vars_long <- vars %>%
    pivot_longer(
      cols = c("SCcat", "SClin", "SW"),
      names_to = "design",
      values_to = "variance"
    ) %>%
    mutate(
      m = ifelse(design=="SW", mSW, mSC),
      psi = ifelse(design=="SW", psiSW, psiSC)
    ) %>%
    select(-c(psiSW, psiSC))
  
  return(vars_long)  
}  

VarSCbasic_SW_multi_line_plot <- function(){

  S1 <- 3
  S2 <- 9

  var_S1 <- Varvals_psi(S1, 1)
  var_S2 <- Varvals_psi(S2, 1)

  allvars <- bind_rows(
    var_S1,
    var_S2
  )

  S.labs <- c("S = 3", "S = 9")
  names(S.labs) <- c(S1, S2)
  
  title <- "Scaled variance of treatment effect estimator"
  p <- ggplot(allvars, aes(x=psi, y=variance, colour=design, linetype=design)) +
    geom_line(size=1.2) +
#    scale_colour_viridis_d() +
    facet_grid(
      . ~ S,
      labeller = labeller(S = S.labs)
    ) +
    xlab(expression(paste("Correlation between cluster-period means, ", psi))) +
    ylab("Scaled variance") +
    labs(title=title, colour="Design",
         linetype="Design") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=16),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
          strip.background = element_rect(
            color="white", fill="white", linetype="solid"
          ),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size=14),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=14), legend.text=element_text(size=14),
          legend.position="bottom")

  ggsave(paste0("plots/var_lines_psi.jpg"),
         p, width=9, height=5, units="in", dpi=600)
  
  return(p)
}

VarSCbasic_SW_multi_line_plot_a <- function(S, K){
  
  rho1 <- 0.01
  rho2 <- 0.1
  rho3 <- 0.25
  m1 <- 10
  m2 <- 100
  
  var_m1_rho1 <- Varvals_a_psi(m1, m1, rho1, S, K)
  var_m2_rho1 <- Varvals_a_psi(m2, m2, rho1, S, K)
  var_m1_rho2 <- Varvals_a_psi(m1, m1, rho2, S, K)
  var_m2_rho2 <- Varvals_a_psi(m2, m2, rho2, S, K)
  var_m1_rho3 <- Varvals_a_psi(m1, m1, rho3, S, K)
  var_m2_rho3 <- Varvals_a_psi(m2, m2, rho3, S, K)
  
  allvars <- bind_rows(
    var_m1_rho1,
    var_m2_rho1,
    var_m1_rho2,
    var_m2_rho2,
    var_m1_rho3,
    var_m2_rho3
  )

  rho.labs <- c(bquote(rho==.(rho1)), bquote(rho==.(rho2)), bquote(rho==.(rho3)))
  names(rho.labs) <- c(rho1, rho2, rho3)
  m.labs <- c(paste("m = ", m1), paste("m = ", m2))
  names(m.labs) <- c(m1, m2)

  title <- "Variance of treatment effect estimator"
  p <- ggplot(allvars, aes(x=r, y=variance, colour=design, linetype=design)) +
    geom_line(size=1.2) +
    facet_grid(
      mSW ~ rho0,
#      labeller = labeller(mSW = m.labs, rho0 = rho.labs)
      labeller = label_bquote(rows=m==.(mSW), cols=rho==.(rho0))
    ) +
#    xlab(expression(paste("Correlation between cluster-period means, ", psi))) +
    xlab(expression(paste("Cluster autocorrelation, ", r))) +
    ylab("Variance") +
    labs(title=title, colour="Design",
         linetype="Design") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=16),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
          strip.background = element_rect(
            color="white", fill="white", linetype="solid"
          ),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size=14),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=14), legend.text=element_text(size=14),
          legend.position="bottom")
  
  ggsave(paste0("plots/variances_lines_S", S, ".jpg"),
         p, width=12, height=7, units="in", dpi=600)
  
  return(p)
}

VarSCbasic_SW_multi_line_plot_a_diffm <- function(S, K){
  
  rho1 <- 0.01
  rho2 <- 0.1
  rho3 <- 0.25

  m1_SW <- 10
  m1_SC <- ((S+1)/2)*m1_SW
  m2_SW <- 100
  m2_SC <- ((S+1)/2)*m2_SW

  var_m1_rho1 <- Varvals_a_psi(m1_SW, m1_SC, rho1, S, K)
  var_m2_rho1 <- Varvals_a_psi(m2_SW, m2_SC, rho1, S, K)
  var_m1_rho2 <- Varvals_a_psi(m1_SW, m1_SC, rho2, S, K)
  var_m2_rho2 <- Varvals_a_psi(m2_SW, m2_SC, rho2, S, K)
  var_m1_rho3 <- Varvals_a_psi(m1_SW, m1_SC, rho3, S, K)
  var_m2_rho3 <- Varvals_a_psi(m2_SW, m2_SC, rho3, S, K)
  
  allvars <- bind_rows(
    var_m1_rho1,
    var_m2_rho1,
    var_m1_rho2,
    var_m2_rho2,
    var_m1_rho3,
    var_m2_rho3
  )
  
  rho.labs <- c(bquote(rho==.(rho1)), bquote(rho==.(rho2)), bquote(rho==.(rho3)))
  names(rho.labs) <- c(rho1, rho2, rho3)
  m.labs <- c(paste("mSC = ", m1_SC, "\nmSW = 10"), paste("mSC = ", m2_SC, "\nmSW = 100"))
  names(m.labs) <- c(m1_SW, m2_SW)

  title <- "Variance of treatment effect estimator"
  p <- ggplot(allvars, aes(x=r, y=variance, colour=design, linetype=design)) +
    geom_line(size=1.2) +
    #    scale_colour_viridis_d() +
    facet_grid(
      mlab ~ rho0,
#      labeller = labeller(mSW = m.labs, rho0 = rho.labs)
      labeller = label_bquote(rows=.(mlab), cols=rho==.(rho0))
    ) +
    #    xlab(expression(paste("Correlation between cluster-period means, ", psi))) +
    xlab(expression(paste("Cluster autocorrelation, ", r))) +
    ylab("Variance") +
    labs(title=title, colour="Design",
         linetype="Design") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=16),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
          strip.background = element_rect(
            color="white", fill="white", linetype="solid"
          ),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size=14),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=14), legend.text=element_text(size=14),
          legend.position="bottom")
  
  ggsave(paste0("plots/var_lines_psi_a_diffm_S", S, ".jpg"),
         p, width=12, height=7, units="in", dpi=600)
  
  return(p)
}

releff_psi <- function(S, K){
  
  psivals <- seq(0, 1, 0.01)
  vars <- data.frame(psi=psivals)
  
  vars <- vars %>%
    mutate(
      SCcat = VarSCcat_no_a(S, K, psi),
      SW = VarSW_no_a(S, K, psi),
      S = S,
      releff = SW/SCcat
    )
  
  return(vars)  
}

Varvals_a_psi <- function(m_SW, m_SC, rho0, S, K){
  
  rvals <- seq(0, 1, 0.01)
  vars <- data.frame(r=rvals)
  
  vars <- vars %>%
    mutate(
      psiSW = psi_corr(m_SW, rho0, r),
      psiSC = psi_corr(m_SC, rho0, r),
      SCcat = VarSCcat(m_SC, S, K, rho0, r),
      SClin = VarSClin(m_SC, S, K, rho0, r),
      SW = VarSW_alt(m_SW, S, K, rho0, r),
      mSW = m_SW,
      mSC = m_SC,
      mlab = paste("mSW = ", mSW, "\nmSC = ", mSC),
      rho0 = rho0,
      S = S
    )
  
  vars_long <- vars %>%
    pivot_longer(
      cols = c("SCcat", "SClin", "SW"),
      names_to = "design",
      values_to = "variance"
    ) %>%
    mutate(
      m = ifelse(design=="SW", mSW, mSC),
      psi = ifelse(design=="SW", psiSW, psiSC)
    ) %>%
    select(-c(psiSW, psiSC))
  
  return(vars_long)  
}  


VarSCbasic_SW_multi_line_plot_psi <- function(){
  
  S1 <- 3
  S2 <- 9
  
  var_S1 <- Varvals_psi(S1, 1)
  var_S2 <- Varvals_psi(S2, 1)
  
  allvars <- bind_rows(
    var_S1,
    var_S2
  )
  
  S.labs <- c("S = 3", "S = 9")
  names(S.labs) <- c(S1, S2)
  
  title <- "Scaled variance of treatment effect estimator"
  p <- ggplot(allvars, aes(x=psi, y=variance, colour=design, linetype=design)) +
    geom_line(size=1.2) +
    #    scale_colour_viridis_d() +
    facet_grid(
      . ~ S,
      labeller = labeller(S = S.labs)
    ) +
    xlab(expression(paste("Correlation between cluster-period means, ", psi))) +
    ylab("Scaled variance") +
    labs(title=title, colour="Design",
         linetype="Design") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=16),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
          strip.background = element_rect(
            color="white", fill="white", linetype="solid"
          ),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size=14),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=14), legend.text=element_text(size=14),
          legend.position="bottom")
  
  ggsave(paste0("plots/var_lines_psi.jpg"),
         p, width=9, height=5, units="in", dpi=600)
  
  return(p)
}


# releffSCSW_grid_multiplot_corr_diffS <- function(S1, S2, K_SW, K_SC, pre_SC, post_SC,
#                                                  corrtype, pereff,
#                                                  fixedscale=FALSE, limits=c(1,2),
#                                                  breaks=seq(1,2,0.5),
#                                                  fullrange=FALSE){
#   # Compare variances of complete SW and staircase designs, for a range of
#   # correlation parameters
#   # Inputs:
#   #  m_SW - number of subjects per cluster-period for SW design
#   #  S_SW - number of unique treatment sequences for SW design
#   #  K_SW - number of times each sequence is repeated for SW design
#   #  m_SC - number of subjects per cluster-period for SC design
#   #  S_SC - number of unique treatment sequences for SC design
#   #  K_SC - number of times each sequence is repeated for SC design
#   #  pre_SC - number of pre-switch measurement periods for SC design
#   #  post_SC - number of post-switch measurement periods for SC design
#   #  corrtype - within-cluster correlation structure type
#   #             (0=block-exchangeable, 1=exponential decay)
#   #  pereff - time period effect type
#   #           ('cat'=categorical period effects, 'lin'=linear period effects)
#   # Output:
#   #  Contour plot of relative variances (vartheta_SC/vartheta_SW)
#   
#   m1 <- 10
#   m2 <- 100
#   relvars_m1_S1 <- gridvals_small(m1, S1, K_SW, m1, S1, K_SC,
#                                   pre_SC, post_SC, corrtype, pereff, fullrange)
#   relvars_m1_S1$m <- m1
#   relvars_m1_S1$S <- S1
#   relvars_m2_S1 <- gridvals_small(m2, S1, K_SW, m2, S1, K_SC,
#                                   pre_SC, post_SC, corrtype, pereff, fullrange)
#   relvars_m2_S1$m <- m2
#   relvars_m2_S1$S <- S1
#   relvars_m1_S2 <- gridvals_small(m1, S2, K_SW, m1, S2, K_SC,
#                                    pre_SC, post_SC, corrtype, pereff, fullrange)
#   relvars_m1_S2$m <- m1
#   relvars_m1_S2$S <- S2
#   relvars_m2_S2 <- gridvals_small(m2, S2, K_SW, m2, S2, K_SC,
#                                    pre_SC, post_SC, corrtype, pereff, fullrange)
#   relvars_m2_S2$m <- m2
#   relvars_m2_S2$S <- S2
#   
#   relvars <- bind_rows(
#     relvars_m1_S1,
#     relvars_m2_S1,
#     relvars_m1_S2,
#     relvars_m2_S2
#   )
#   
#   if(fixedscale==TRUE){
#     fillopt <- scale_fill_continuous_divergingx(name="Relative efficiency ",
#                                                 palette='RdBu', mid=1, rev=FALSE,
#                                                 limits=limits, breaks=breaks) #palette='RdBu' 'Temps' 'RdYlBu'
#   }else{
#     fillopt <- scale_fill_continuous_divergingx(name="Relative efficiency ",
#                                                 palette='RdBu', mid=1, rev=FALSE)
#   }
#   
#   S.labs <- c(paste("S = ", S1), paste("S = ", S2))
#   names(S.labs) <- c(S1, S2)
#   m.labs <- c(paste("m = ", m1), paste("m = ", m2))
#   names(m.labs) <- c(m1, m2)
#   
#   KSW <- ifelse(K_SW==1, "", K_SW)
#   KSC <- ifelse(K_SC==1, "", K_SC)
#   title <- bquote(paste("Relative efficiency, ",
#                         var(hat(theta))[paste("SW(S,", .(KSW), "k)")]/
#                           var(hat(theta))[paste("SC(S,", .(KSC), "k,", .(pre_SC), ",", .(post_SC), ")")]))
# 
#   p <- ggplot(relvars, aes(x=r, y=rho)) +
#     geom_tile(aes(fill=releffSCSW)) +
#     fillopt +
#     facet_grid(
#       m ~ S,
#       labeller = labeller(m = m.labs, S = S.labs)
#     ) +
#     scale_x_continuous(expand=c(0,0)) +
#     scale_y_continuous(expand=c(0,0)) +
#     theme(aspect.ratio=3/8,
#           legend.key.width = unit(1.5, "cm"),
#           legend.title=element_text(size=14), legend.text=element_text(size=14),
#           legend.position="bottom",
#           plot.title=element_text(hjust=0.5, size=16),
#           axis.title=element_text(size=14), axis.text=element_text(size=14),
#           strip.background = element_rect(
#             color="white", fill="white", linetype="solid"
#           ),
#           strip.text.x = element_text(size = 12),
#           strip.text.y = element_text(size=12)) +
#     coord_fixed() + xlab(expression(paste("Cluster autocorrelation, ", r))) + ylab(expression(paste("Within-period ICC, ", rho))) +
#     ggtitle(title)
# 
#   rng <- ifelse(fullrange, "full", "restricted")
#   corrstruct <- ifelse(corrtype==0, "BE", "DTD")
#   
#   ggsave(paste0("plots/releff_SC_S1", pre_SC, post_SC, "_vs_SW_S1_",
#                 corrstruct, "_", pereff, "_", rng, "_S_", S1, "_", S2, ".jpg"),
#          p, width=9, height=5, units="in", dpi=800)
#   return(p)
# }

releffSCSW_grid_multiplot_corr_diffS <- function(S1, S2, K_SW, K_SC,
                                                 pre_SC, post_SC,
                                                 corrtype, pereff,
                                                 fixedscale=FALSE, limits=c(1,2),
                                                 breaks=seq(1,2,0.5),
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
  relvars_m1_S1 <- gridvals_small(m1, S1, K_SW, m1, S1, K_SC,
                                  pre_SC, post_SC, corrtype, pereff, fullrange)
  relvars_m1_S1$m <- m1
  relvars_m1_S1$S <- S1
  relvars_m2_S1 <- gridvals_small(m2, S1, K_SW, m2, S1, K_SC,
                                  pre_SC, post_SC, corrtype, pereff, fullrange)
  relvars_m2_S1$m <- m2
  relvars_m2_S1$S <- S1
  relvars_m1_S2 <- gridvals_small(m1, S2, K_SW, m1, S2, K_SC,
                                   pre_SC, post_SC, corrtype, pereff, fullrange)
  relvars_m1_S2$m <- m1
  relvars_m1_S2$S <- S2
  relvars_m2_S2 <- gridvals_small(m2, S2, K_SW, m2, S2, K_SC,
                                   pre_SC, post_SC, corrtype, pereff, fullrange)
  relvars_m2_S2$m <- m2
  relvars_m2_S2$S <- S2
  
  relvars <- bind_rows(
    relvars_m1_S1,
    relvars_m2_S1,
    relvars_m1_S2,
    relvars_m2_S2
  )
  
  if(fixedscale==TRUE){
    fillopt <- scale_fill_continuous_divergingx(name="Relative efficiency ",
                                                palette='RdBu', mid=1, rev=FALSE,
                                                limits=limits, breaks=breaks) #palette='RdBu' 'Temps' 'RdYlBu'
  }else{
    fillopt <- scale_fill_continuous_divergingx(name="Relative efficiency ",
                                                palette='RdBu', mid=1, rev=FALSE)
  }
  
  S.labs <- c(paste("S = ", S1), paste("S = ", S2))
  names(S.labs) <- c(S1, S2)
  m.labs <- c(paste("m = ", m1), paste("m = ", m2))
  names(m.labs) <- c(m1, m2)
  
  KSW <- ifelse(K_SW==1, "", K_SW)
  KSC <- ifelse(K_SC==1, "", K_SC)
  title <- bquote(paste("Relative efficiency, ",
                        var(hat(theta))[paste("SW(S,", .(KSW), "k)")]/
                          var(hat(theta))[paste("SC(S,", .(KSC), "k,", .(pre_SC), ",", .(post_SC), ")")]))

  p <- ggplot(relvars, aes(x=r, y=rho)) +
    geom_tile(aes(fill=releffSCSW)) +
    fillopt +
    facet_grid(
      m ~ S,
      labeller = labeller(m = m.labs, S = S.labs)
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
    ggtitle(title)

  rng <- ifelse(fullrange, "full", "restricted")
  corrstruct <- ifelse(corrtype==0, "BE", "DTD")
  
  ggsave(paste0("plots/releff_SC_S1", pre_SC, post_SC, "_vs_SW_S1_",
                corrstruct, "_", pereff, "_", rng, "_S_", S1, "_", S2, ".jpg"),
         p, width=9, height=5, units="in", dpi=800)
  ggsave(paste0("plots/releff_SC_S1", pre_SC, post_SC, "_vs_SW_S1_",
                corrstruct, "_", pereff, "_", rng, "_S_", S1, "_", S2, ".pdf"),
         p, width=9, height=5, units="in", dpi=600)
  
  return(p)
}

releffSCSW_grid_multiplot_diffm_diffS <- function(S1, S2, pre_SC, post_SC,
                                                  corrtype, pereff,
                                                  fixedscale=FALSE, limits=c(1,2),
                                                  breaks=seq(1,2,0.5),
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
  m1_S1_SC <- ((S1+1)/2)*m1_SW
  m1_S2_SC <- ((S2+1)/2)*m1_SW
  m2_SW <- 100
  m2_S1_SC <- ((S1+1)/2)*m2_SW
  m2_S2_SC <- ((S2+1)/2)*m2_SW
  relvars_m1_S1 <- gridvals_small(m1_SW, S1, 1, m1_S1_SC, S1, 1,
                                  pre_SC, post_SC, corrtype, pereff, fullrange)
  relvars_m1_S1$m <- m1_SW
  relvars_m1_S1$S <- S1
  relvars_m2_S1 <- gridvals_small(m2_SW, S1, 1, m2_S1_SC, S1, 1,
                                  pre_SC, post_SC, corrtype, pereff, fullrange)
  relvars_m2_S1$m <- m2_SW
  relvars_m2_S1$S <- S1
  relvars_m1_S2 <- gridvals_small(m1_SW, S2, 1, m1_S2_SC, S2, 1,
                                   pre_SC, post_SC, corrtype, pereff, fullrange)
  relvars_m1_S2$m <- m1_SW
  relvars_m1_S2$S <- S2
  relvars_m2_S2 <- gridvals_small(m2_SW, S2, 1, m2_S2_SC, S2, 1,
                                   pre_SC, post_SC, corrtype, pereff, fullrange)
  relvars_m2_S2$m <- m2_SW
  relvars_m2_S2$S <- S2
  relvars <- bind_rows(
    relvars_m1_S1,
    relvars_m2_S1,
    relvars_m1_S2,
    relvars_m2_S2
  )
  
  if(fixedscale==TRUE){
    fillopt <- scale_fill_continuous_divergingx(name="Relative efficiency ",
                                                palette='RdBu', mid=1, rev=FALSE,
                                                limits=limits, breaks=breaks)
  }else{
    fillopt <- scale_fill_continuous_divergingx(name="Relative efficiency ",
                                                palette='RdBu', mid=1, rev=FALSE)
  }
  
  S.labs <- c(paste("S = ", S1), paste("S = ", S2))
  names(S.labs) <- c(S1, S2)
  m.labs <- c("mSW=10 \nmSC=5*(S+1)", "mSW=100 \nmSC=50*(S+1)")
  names(m.labs) <- c(m1_SW, m2_SW)
  
  p <- ggplot(relvars, aes(x=r, y=rho)) + 
    geom_tile(aes(fill=releffSCSW)) +
    facet_grid(
      m ~ S,
      labeller = labeller(m = m.labs, S = S.labs)
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
          strip.text.x = element_text(size=10),
          strip.text.y = element_text(size=10)) +
    coord_fixed() + xlab(expression(paste("Cluster autocorrelation, ", r))) + ylab(expression(paste("Within-period ICC, ", rho))) +
    ggtitle(bquote(paste("Relative efficiency, ",
                         var(hat(theta))["SW(S,k)"]/
                           var(hat(theta))[paste("SC(S,k,", .(pre_SC), ",", .(post_SC), ")")])))
  
  rng <- ifelse(fullrange, "full", "restricted")
  corrstruct <- ifelse(corrtype==0, "BE", "DTD")
  
  ggsave(paste0("plots/releff_diffm_SC_S1", pre_SC, post_SC, "_vs_SW_S1",
                corrstruct, "_", pereff, "_", rng, "_S_", S1, "_", S2, ".jpg"),
         p, width=9, height=5, units="in", dpi=800)
  ggsave(paste0("plots/releff_diffm_SC_S1", pre_SC, post_SC, "_vs_SW_S1",
                corrstruct, "_", pereff, "_", rng, "_S_", S1, "_", S2, ".pdf"),
         p, width=9, height=5, units="in", dpi=600)
  
  return(p)
}

releff_SW_extendedSC_diffS <- function(S1, S2, K_SW,
                                       pre_SC, post_SC,
                                       corrtype, pereff,
                                       fixedscale=FALSE, limits=c(1,2),
                                       breaks=seq(1,2,0.5),
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
  K_S1_SC <- K_SW*(S1+1)/(pre_SC + post_SC)
  K_S2_SC <- K_SW*(S2+1)/(pre_SC + post_SC)
  relvars_m1_S1 <- gridvals_small(m1, S1, K_SW, m1, S1, K_S1_SC,
                                  pre_SC, post_SC, corrtype, pereff, fullrange)
  relvars_m1_S1$m <- m1
  relvars_m1_S1$S <- S1
  relvars_m2_S1 <- gridvals_small(m2, S1, K_SW, m2, S1, K_S1_SC,
                                  pre_SC, post_SC, corrtype, pereff, fullrange)
  relvars_m2_S1$m <- m2
  relvars_m2_S1$S <- S1
  relvars_m1_S2 <- gridvals_small(m1, S2, K_SW, m1, S2, K_S2_SC,
                                  pre_SC, post_SC, corrtype, pereff, fullrange)
  relvars_m1_S2$m <- m1
  relvars_m1_S2$S <- S2
  relvars_m2_S2 <- gridvals_small(m2, S2, K_SW, m2, S2, K_S2_SC,
                                  pre_SC, post_SC, corrtype, pereff, fullrange)
  relvars_m2_S2$m <- m2
  relvars_m2_S2$S <- S2
  
  relvars <- bind_rows(
    relvars_m1_S1,
    relvars_m2_S1,
    relvars_m1_S2,
    relvars_m2_S2
  )
  
  if(fixedscale==TRUE){
    fillopt <- scale_fill_continuous_divergingx(name="Relative efficiency ",
                                                palette='RdBu', mid=1, rev=FALSE,
                                                limits=limits, breaks=breaks) #palette='RdBu' 'Temps' 'RdYlBu'
  }else{
    fillopt <- scale_fill_continuous_divergingx(name="Relative efficiency ",
                                                palette='RdBu', mid=1, rev=FALSE)
  }
  
  S.labs <- c(paste("S = ", S1, ", q = ", K_S1_SC), paste("S = ", S2, ", q = ", K_S2_SC))
  names(S.labs) <- c(S1, S2)
  m.labs <- c(paste("m = ", m1), paste("m = ", m2))
  names(m.labs) <- c(m1, m2)
  
  KSW <- ifelse(K_SW==1, "", K_SW)
  title <- bquote(paste("Relative efficiency, ",
                        var(hat(theta))[paste("SW(S,", .(KSW), "k)")]/
                          var(hat(theta))[paste("SC(S,qk,", .(pre_SC), ",", .(post_SC), ")")]))
  
  p <- ggplot(relvars, aes(x=r, y=rho)) +
    geom_tile(aes(fill=releffSCSW)) +
    fillopt +
    facet_grid(
      m ~ S,
      labeller = labeller(m = m.labs, S = S.labs)
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
    ggtitle(title)
  
  rng <- ifelse(fullrange, "full", "restricted")
  corrstruct <- ifelse(corrtype==0, "BE", "DTD")
  
  ggsave(paste0("plots/releff_extendedSC_Sq", pre_SC, post_SC, "_vs_SW_S", K_SW,
                corrstruct, "_", pereff, "_", rng, "_S_", S1, "_", S2, ".jpg"),
         p, width=9, height=5, units="in", dpi=800)
  ggsave(paste0("plots/releff_extendedSC_Sq", pre_SC, post_SC, "_vs_SW_S", K_SW,
                corrstruct, "_", pereff, "_", rng, "_S_", S1, "_", S2, ".pdf"),
         p, width=9, height=5, units="in", dpi=600)
  
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
