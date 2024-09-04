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
library(metR)

gridvals_small <- function(m_SW, S_SW, reps_SW, m_SC, S_SC, reps_SC,
                           corrtype, pereff, fullrange=FALSE, rhou=0, imp=FALSE){
  # Compare variances of complete SW and staircase designs, for a range of
  # correlation parameters
  # Inputs:
  #  m_SW - number of subjects per cluster-period for SW design
  #  S_SW - number of unique treatment sequences for SW design
  #  reps_SW - number of times each sequence is repeated for SW design
  #  m_SC - number of subjects per cluster-period for SC design
  #  S_SC - number of unique treatment sequences for SC design
  #  reps_SC - number of times each sequence is repeated for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  #  fullrange - Use full correlation parameter ranges (T/F)
  #  rhou - correlation between a participant's repeated measurements
  #          (default of 0 means each participant is measured just once)
  #  imp - Include implementation periods in design (T/F)
  # Output:
  #  Relative efficiencies (vartheta_SW/vartheta_SC) for a range of rho and r
  #  values, for a particular value of rhou
  
  if(fullrange==TRUE){
    rhovals <- seq(0.01, 0.99, 0.01)
    rvals <- seq(0.05, 1, 0.05)
  }else{
    rhovals <- seq(0.01, 0.25, 0.01)
    rvals <- seq(0.2, 1, 0.05)
  }
  vars <- expand.grid(rho=rhovals, r=rvals)
  
  vars$varSW <- with(
    vars,
    sapply(1:nrow(vars), function(j){
      CRTvartheta(m_SW, SWdesmat(S_SW, reps_SW, imp),
               rho[j], r[j], corrtype, pereff, rhou)
      }
    )
  )
  
  if(imp){ # If including implementation periods, must use general variance expression
    vars$varSC <- with(
      vars,
      sapply(1:nrow(vars), function(j){
        CRTvartheta(m_SC, SCdesmat(S_SC, reps_SC, 1, 1, imp),
                    rho[j], r[j], corrtype, pereff, rhou)
        }
      )
    )
  }else{ # If not including implementation periods, can use explicit variance expressions
    if(pereff=="cat"){
      vars <- vars %>%
        mutate(
          varSC = VarSCcat(m_SC, S_SC, reps_SC, rho, r, rhou),
        )
    }else if(pereff=="lin"){
      vars <- vars %>%
        mutate(
          varSC = VarSClin(m_SC, S_SC, reps_SC, rho, r, rhou),
        )
    }
  }
  
  vars <- vars %>%
    mutate(
      releffSCSW = varSW/varSC
    )
  return(vars)
}

corr_grid_multiplot <- function(fixedscale=TRUE, limits=c(0,1),
                                breaks=seq(0,1,0.2), fullrange=FALSE, rhou=0){
  
  if(fullrange==TRUE){
    rhovals <- seq(0.01, 0.99, 0.01)
    rvals <- seq(0.05, 1, 0.05)
    mvals <- c(10, 100)
  }else{
    rhovals <- seq(0.01, 0.25, 0.01)
    rvals <- seq(0.2, 1, 0.05)
    mvals <- c(10, 100)
  }
  vars <- expand.grid(rho=rhovals, r=rvals, mvals=mvals, rhou=rhou)
  
  vars$m <- ifelse(vars$mvals==10, "m = 10", "m = 100")
  
  vars$corr <- with(
    vars,
    sapply(1:nrow(vars), function(j){
      psi_corr(mvals[j], rho[j], r[j], rhou[j])
    }
    )
  )
  
  if(fixedscale==TRUE){
    fillopt <- scale_fill_distiller(name="Correlation", 
                                    palette="YlOrRd", direction=1,
                                    limits=limits, breaks=breaks)
  }else{
    fillopt <- scale_fill_distiller(name="Correlation",
                                    palette="YlOrRd", direction=1)
  }
  
  if(rhou!=0){
    stitle <- bquote(paste("Closed cohort with ", rho[u], "=", .(rhou)))
  }else{
    stitle <- NULL
  }

  p <- ggplot(vars, aes(x=r, y=rho, z=corr)) + 
    geom_tile(aes(fill=corr)) +
    geom_contour2(color="black", breaks=breaks) +
    geom_text_contour(stroke = 0.15, size=3, rotate=FALSE, breaks=breaks,
                      skip=0, label.placer=label_placer_fraction(0.5)) +
    facet_grid(m ~ .) +
    fillopt +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(aspect.ratio=3/8,
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=14), legend.text=element_text(size=14),
          legend.position="bottom",
          plot.title=element_text(hjust=0.5, size=16),
          plot.subtitle=element_text(hjust=0.5, size=14),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
          strip.background = element_rect(
            color="white", fill="white", linetype="solid"
          ),
          strip.text.x = element_text(size = 10),
          strip.text.y = element_text(size=10)) +
    coord_fixed() + xlab(expression(paste("Cluster autocorrelation, ", r))) + ylab(expression(paste("Within-period ICC, ", rho))) +
    ggtitle(expression(paste("Correlation between cluster-period means, ", psi)),
            subtitle=stitle)
  
  rng <- ifelse(fullrange, "full", "restricted")
  rhouval <- ifelse((0 < rhou && rhou < 1), strsplit(as.character(rhou), "\\.")[[1]][2],
                    as.character(rhou))
  
  ggsave(paste0("plots/corr_m10_m100_", rng, "_rhou_", rhouval, ".jpg"),
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
                                              limits=c(0,1), breaks=seq(0,1,0.25))

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

releffSCSW_grid_multiplot_corr_diffS <- function(S1, S2, K_SW, K_SC,
                                                 corrtype, pereff,
                                                 fixedscale=FALSE, limits=c(0,1),
                                                 breaks=seq(0,1,0.2),
                                                 fullrange=FALSE, rhou=0, imp=FALSE){
  # Compare variances of complete SW designs and embedded basic staircase designs,
  # for a range of correlation parameters, and for different numbers of sequences
  # and cluster-period sizes
  # Inputs:
  #  S1 - number of unique treatment sequences for left column
  #  S2 - number of unique treatment sequences for right column
  #  K_SW - number of times each sequence is repeated for SW design
  #  K_SC - number of times each sequence is repeated for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  #  fixedscale - Use particular limits and breaks for colour scale (T/F)
  #  limits - Range of values for colour scale shading
  #  breaks - Indices for colour scale legend and contour lines
  #  fullrange - Use full correlation parameter ranges (T/F)
  #  rhou - correlation between a participant's repeated measurements
  #          (default of 0 means each participant is measured just once)
  #  imp - Include implementation periods in design (T/F)
  # Output:
  #  Contour plot of relative efficiencies (vartheta_SW/vartheta_SC)
  
  m1 <- 10
  m2 <- 100
  relvars_m1_S1 <- gridvals_small(m1, S1, K_SW, m1, S1, K_SC,
                                  corrtype, pereff, fullrange, rhou, imp)
  relvars_m1_S1$m <- m1
  relvars_m1_S1$S <- S1
  relvars_m2_S1 <- gridvals_small(m2, S1, K_SW, m2, S1, K_SC,
                                  corrtype, pereff, fullrange, rhou, imp)
  relvars_m2_S1$m <- m2
  relvars_m2_S1$S <- S1
  relvars_m1_S2 <- gridvals_small(m1, S2, K_SW, m1, S2, K_SC,
                                  corrtype, pereff, fullrange, rhou, imp)
  relvars_m1_S2$m <- m1
  relvars_m1_S2$S <- S2
  relvars_m2_S2 <- gridvals_small(m2, S2, K_SW, m2, S2, K_SC,
                                  corrtype, pereff, fullrange, rhou, imp)
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
                                                limits=limits, breaks=breaks)
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
                        var(hat(theta))[paste("SW(S,", .(KSW), "K)")]/
                          var(hat(theta))[paste("SC(S,", .(KSC), "K,1,1)")]))
  if(rhou==0 && imp){
    stitle <- "Including implementation periods"
  }else if(rhou!=0 && imp){
    stitle <- bquote(paste("Closed cohort with ", rho[u], "=", .(rhou),
                           ", including implementation periods"))
  }else if(rhou!=0 && isFALSE(imp)){
    stitle <- bquote(paste("Closed cohort with ", rho[u], "=", .(rhou)))
  }else{
    stitle <- NULL
  }

  p <- ggplot(relvars, aes(x=r, y=rho, z=releffSCSW)) +
    geom_tile(aes(fill=releffSCSW)) +
#    geom_contour_fill(breaks=breaks) +
    geom_contour2(color="black", breaks=breaks) +
    geom_text_contour(stroke = 0.15, size=3, rotate=FALSE, breaks=breaks,
                      skip=0, label.placer=label_placer_fraction(0.5)) +
    fillopt +
    facet_grid(
      m ~ S,
      labeller = labeller(m = m.labs, S = S.labs)
    ) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(aspect.ratio=3/8,
          panel.spacing = unit(1.5, "lines"),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=14), legend.text=element_text(size=14),
          legend.position="bottom",
          plot.title=element_text(hjust=0.5, size=16),
          plot.subtitle=element_text(hjust=0.5, size=14),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
          strip.background = element_rect(
            color="white", fill="white", linetype="solid"
          ),
          strip.text.x = element_text(size=12),
          strip.text.y = element_text(size=12)) +
    coord_fixed() + xlab(expression(paste("Cluster autocorrelation, ", r))) +
    ylab(expression(paste("Within-period ICC, ", rho))) +
    ggtitle(title, subtitle=stitle)

  rng <- ifelse(fullrange, "full", "restricted")
  corrstruct <- ifelse(corrtype==0, "BE", "DTD")
  rhouval <- ifelse((0 < rhou && rhou < 1), strsplit(as.character(rhou), "\\.")[[1]][2],
                    as.character(rhou))
  impper <- ifelse(imp, "_implementation", "")
  ggsave(paste0("plots/releff_SC_vs_SW_",
                corrstruct, "_", pereff, "_", rng, "_S_", S1, "_", S2, "_rhou_",
                rhouval, impper, ".jpg"), p, width=9, height=5, units="in", dpi=800)
  ggsave(paste0("plots/releff_SC_vs_SW_",
                corrstruct, "_", pereff, "_", rng, "_S_", S1, "_", S2, "_rhou_",
                rhouval, impper, ".pdf"), p, width=9, height=5, units="in", dpi=600)
  return(p)
}

releffSCSW_grid_multiplot_diffm_diffS <- function(S1, S2,
                                                  corrtype, pereff,
                                                  fixedscale=FALSE, limits=c(0,2),
                                                  breaks=seq(0,2,0.2),
                                                  fullrange=FALSE, rhou=0, imp=FALSE){
  # Compare variances of complete SW designs and basic staircase designs with
  # larger cluster-period sizes (the same total number of participants), for a
  # range of correlation parameters, and for different numbers of sequences and
  # cluster-period sizes
  # Inputs:
  #  S1 - number of unique treatment sequences for left column
  #  S2 - number of unique treatment sequences for right column
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  #  fixedscale - Use particular limits and breaks for colour scale (T/F)
  #  limits - Range of values for colour scale shading
  #  breaks - Indices for colour scale legend and contour lines
  #  fullrange - Use full correlation parameter ranges (T/F)
  #  rhou - correlation between a participant's repeated measurements
  #          (default of 0 means each participant is measured just once)
  # Output:
  #  Contour plot of relative efficiencies (vartheta_SW/vartheta_SC)

  m1_SW <- 10
  m1_S1_SC <- ((S1+1)/2)*m1_SW
  m1_S2_SC <- ((S2+1)/2)*m1_SW
  m2_SW <- 100
  m2_S1_SC <- ((S1+1)/2)*m2_SW
  m2_S2_SC <- ((S2+1)/2)*m2_SW
  relvars_m1_S1 <- gridvals_small(m1_SW, S1, 1, m1_S1_SC, S1, 1,
                                  corrtype, pereff, fullrange, rhou, imp)
  relvars_m1_S1$m <- m1_SW
  relvars_m1_S1$S <- S1
  relvars_m2_S1 <- gridvals_small(m2_SW, S1, 1, m2_S1_SC, S1, 1,
                                  corrtype, pereff, fullrange, rhou, imp)
  relvars_m2_S1$m <- m2_SW
  relvars_m2_S1$S <- S1
  relvars_m1_S2 <- gridvals_small(m1_SW, S2, 1, m1_S2_SC, S2, 1,
                                  corrtype, pereff, fullrange, rhou, imp)
  relvars_m1_S2$m <- m1_SW
  relvars_m1_S2$S <- S2
  relvars_m2_S2 <- gridvals_small(m2_SW, S2, 1, m2_S2_SC, S2, 1,
                                  corrtype, pereff, fullrange, rhou, imp)
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
  
  title <- bquote(paste("Relative efficiency, ",
                        var(hat(theta))["SW(S,K)"]/
                          var(hat(theta))["SC(S,K,1,1)"]))
  if(rhou==0 && imp){
    stitle <- "Including implementation periods"
  }else if(rhou!=0 && imp){
    stitle <- bquote(paste("Closed cohort with ", rho[u], "=", .(rhou),
                           ", including implementation periods"))
  }else if(rhou!=0 && isFALSE(imp)){
    stitle <- bquote(paste("Closed cohort with ", rho[u], "=", .(rhou)))
  }else{
    stitle <- NULL
  }
  
  p <- ggplot(relvars, aes(x=r, y=rho, z=releffSCSW)) + 
    geom_tile(aes(fill=releffSCSW)) +
#    geom_contour_fill(breaks=breaks) +
    geom_contour2(color="black", breaks=breaks) +
    geom_text_contour(stroke=0.15, size=3, rotate=FALSE, breaks=breaks, skip=0,
                      label.placer=label_placer_fraction(0.5)) +
    facet_grid(
      m ~ S,
      labeller = labeller(m = m.labs, S = S.labs)
    ) +
    fillopt +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(aspect.ratio=3/8,
          panel.spacing = unit(1.5, "lines"),
          legend.key.width = unit(2.5, "cm"),
          legend.title=element_text(size=14), legend.text=element_text(size=14),
          legend.position="bottom",
          plot.title=element_text(hjust=0.5, size=16),
          plot.subtitle=element_text(hjust=0.5, size=14),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
          strip.background = element_rect(
            color="white", fill="white", linetype="solid"
          ),
          strip.text.x = element_text(size=12),
          strip.text.y = element_text(size=12)) +
    coord_fixed() + xlab(expression(paste("Cluster autocorrelation, ", r))) +
    ylab(expression(paste("Within-period ICC, ", rho))) +
    ggtitle(title, subtitle=stitle)
  
  rng <- ifelse(fullrange, "full", "restricted")
  corrstruct <- ifelse(corrtype==0, "BE", "DTD")
  rhouval <- ifelse((0 < rhou && rhou < 1), strsplit(as.character(rhou), "\\.")[[1]][2],
                    as.character(rhou))
  impper <- ifelse(imp, "_implementation", "")
  ggsave(paste0("plots/releff_diffm_SC_vs_SW_",
                corrstruct, "_", pereff, "_", rng, "_S_", S1, "_", S2, "_rhou_",
                rhouval, impper, ".jpg"), p, width=9, height=5, units="in", dpi=800)
  ggsave(paste0("plots/releff_diffm_SC_vs_SW_",
                corrstruct, "_", pereff, "_", rng, "_S_", S1, "_", S2, "_rhou_",
                rhouval, impper, ".pdf"), p, width=9, height=5, units="in", dpi=600)
  return(p)
}

releff_SW_extendedSC_diffS <- function(S1, S2, K_SW,
                                       corrtype, pereff,
                                       fixedscale=FALSE, limits=c(0,2.5),
                                       breaks=seq(0,2.5,0.5),
                                       fullrange=FALSE, rhou=0, imp=FALSE){
  # Compare variances of complete SW designs and basic staircase designs with
  # more clusters (the same total number of participants), for a range of
  # correlation parameters, and for different numbers of sequences and
  # cluster-period sizes
  # Inputs:
  #  S1 - number of unique treatment sequences for left column
  #  S2 - number of unique treatment sequences for right column
  #  K_SW - number of times each sequence is repeated for SW design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  #  fixedscale - Use particular limits and breaks for colour scale (T/F)
  #  limits - Range of values for colour scale shading
  #  breaks - Indices for colour scale legend and contour lines
  #  fullrange - Use full correlation parameter ranges (T/F)
  #  rhou - correlation between a participant's repeated measurements
  #          (default of 0 means each participant is measured just once)
  # Output:
  #  Contour plot of relative efficiencies (vartheta_SW/vartheta_SC)

  m1 <- 10
  m2 <- 100
  K_S1_SC <- K_SW*(S1+1)/2
  K_S2_SC <- K_SW*(S2+1)/2
  relvars_m1_S1 <- gridvals_small(m1, S1, K_SW, m1, S1, K_S1_SC,
                                  corrtype, pereff, fullrange, rhou, imp)
  relvars_m1_S1$m <- m1
  relvars_m1_S1$S <- S1
  relvars_m2_S1 <- gridvals_small(m2, S1, K_SW, m2, S1, K_S1_SC,
                                  corrtype, pereff, fullrange, rhou, imp)
  relvars_m2_S1$m <- m2
  relvars_m2_S1$S <- S1
  relvars_m1_S2 <- gridvals_small(m1, S2, K_SW, m1, S2, K_S2_SC,
                                  corrtype, pereff, fullrange, rhou, imp)
  relvars_m1_S2$m <- m1
  relvars_m1_S2$S <- S2
  relvars_m2_S2 <- gridvals_small(m2, S2, K_SW, m2, S2, K_S2_SC,
                                  corrtype, pereff, fullrange, rhou, imp)
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
                                                limits=limits, breaks=breaks)
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
                        var(hat(theta))[paste("SW(S,", .(KSW), "K)")]/
                          var(hat(theta))["SC(S,qK,1,1)"]))
  if(rhou==0 && imp){
    stitle <- "Including implementation periods"
  }else if(rhou!=0 && imp){
    stitle <- bquote(paste("Closed cohort with ", rho[u], "=", .(rhou),
                           ", including implementation periods"))
  }else if(rhou!=0 && isFALSE(imp)){
    stitle <- bquote(paste("Closed cohort with ", rho[u], "=", .(rhou)))
  }else{
    stitle <- NULL
  }
  
  p <- ggplot(relvars, aes(x=r, y=rho, z=releffSCSW)) +
    geom_tile(aes(fill=releffSCSW)) +
#    geom_contour_fill(breaks=breaks) +
    geom_contour2(color="black", breaks=breaks) +
    geom_text_contour(stroke=0.15, size=3, rotate=FALSE, breaks=breaks, skip=0,
                      label.placer=label_placer_fraction(0.5)) +
    fillopt +
    facet_grid(
      m ~ S,
      labeller = labeller(m = m.labs, S = S.labs)
    ) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(aspect.ratio=3/8,
          panel.spacing = unit(1.5, "lines"),
          legend.key.width = unit(2.5, "cm"),
          legend.title=element_text(size=14), legend.text=element_text(size=14),
          legend.position="bottom",
          plot.title=element_text(hjust=0.5, size=16),
          plot.subtitle=element_text(hjust=0.5, size=14),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
          strip.background = element_rect(
            color="white", fill="white", linetype="solid"
          ),
          strip.text.x = element_text(size=12),
          strip.text.y = element_text(size=12)) +
    coord_fixed() + xlab(expression(paste("Cluster autocorrelation, ", r))) +
    ylab(expression(paste("Within-period ICC, ", rho))) +
    ggtitle(title, subtitle=stitle)
  
  rng <- ifelse(fullrange, "full", "restricted")
  corrstruct <- ifelse(corrtype==0, "BE", "DTD")
  rhouval <- ifelse((0 < rhou && rhou < 1), strsplit(as.character(rhou), "\\.")[[1]][2],
                    as.character(rhou))
  impper <- ifelse(imp, "_implementation", "")
  ggsave(paste0("plots/releff_extendedSC_Sq11_vs_SW_S", K_SW, "_",
                corrstruct, "_", pereff, "_", rng, "_S_", S1, "_", S2, "_rhou_",
                rhouval, impper, ".jpg"), p, width=9, height=5, units="in", dpi=800)
  ggsave(paste0("plots/releff_extendedSC_Sq11_vs_SW_S", K_SW, "_",
                corrstruct, "_", pereff, "_", rng, "_S_", S1, "_", S2, "_rhou_",
                rhouval, impper, ".pdf"), p, width=9, height=5, units="in", dpi=600)
  return(p)
}

releffSCSW_grid_multiplot_corr_singleS <- function(S, K, corrtype, pereff,
                                                   fixedscale=FALSE, limits=c(0,1),
                                                   breaks=seq(0,1,0.2),
                                                   fullrange=FALSE, rhou=0, imp=FALSE){
  # Compare variances of complete SW designs and embedded basic staircase designs,
  # for a range of correlation parameters, and for different cluster-period sizes
  # Inputs:
  #  S - number of unique treatment sequences
  #  K - number of clusters randomised to each sequence
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  #  fixedscale - Use particular limits and breaks for colour scale (T/F)
  #  limits - Range of values for colour scale shading
  #  breaks - Indices for colour scale legend and contour lines
  #  fullrange - Use full correlation parameter ranges (T/F)
  #  rhou - correlation between a participant's repeated measurements
  #          (default of 0 means each participant is measured just once)
  #  imp - Include implementation periods in design (T/F)
  # Output:
  #  Contour plot of relative efficiencies (vartheta_SW/vartheta_SC)
  
  m1 <- 10
  m2 <- 100
  relvars_m1 <- gridvals_small(m1, S, K, m1, S, K,
                                  corrtype, pereff, fullrange, rhou, imp)
  relvars_m1$m <- m1
  relvars_m1$S <- S
  relvars_m2 <- gridvals_small(m2, S, K, m2, S, K,
                                  corrtype, pereff, fullrange, rhou, imp)
  relvars_m2$m <- m2
  relvars_m2$S <- S

  relvars <- bind_rows(
    relvars_m1,
    relvars_m2
  )
  
  if(fixedscale==TRUE){
    fillopt <- scale_fill_continuous_divergingx(name="Relative efficiency ",
                                                palette='RdBu', mid=1, rev=FALSE,
                                                limits=limits, breaks=breaks)
  }else{
    fillopt <- scale_fill_continuous_divergingx(name="Relative efficiency ",
                                                palette='RdBu', mid=1, rev=FALSE)
  }
  
  m.labs <- c(paste("m = ", m1), paste("m = ", m2))
  names(m.labs) <- c(m1, m2)
  S.labs <- paste("S = ", S)
  names(S.labs) <- S
  
  Kval <- ifelse(K==1, "", K)
  title <- bquote(paste("Relative efficiency, ",
                        var(hat(theta))[paste("SW(", .(S), ",", .(Kval), "K)")]/
                          var(hat(theta))[paste("SC(", .(S), ",", .(Kval), "K,1,1)")]))
  if(rhou==0 && imp){
    stitle <- "Including implementation periods"
  }else if(rhou!=0 && imp){
    stitle <- bquote(paste("Closed cohort with ", rho[u], "=", .(rhou),
                           ", including implementation periods"))
  }else if(rhou!=0 && isFALSE(imp)){
    stitle <- bquote(paste("Closed cohort with ", rho[u], "=", .(rhou)))
  }else{
    stitle <- NULL
  }
  
  p <- ggplot(relvars, aes(x=r, y=rho, z=releffSCSW)) +
    geom_tile(aes(fill=releffSCSW)) +
    #    geom_contour_fill(breaks=breaks) +
    geom_contour2(color="black", breaks=breaks) +
    geom_text_contour(stroke = 0.15, size=3, rotate=FALSE, breaks=breaks,
                      skip=0, label.placer=label_placer_fraction(0.5)) +
    fillopt +
    facet_grid(
      m ~ S,
      labeller = labeller(m = m.labs, S = S.labs)
    ) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(aspect.ratio=3/8,
#          panel.spacing = unit(1.5, "lines"),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=14), legend.text=element_text(size=14),
          legend.position="bottom",
          plot.title=element_text(hjust=0.5, size=16),
          plot.subtitle=element_text(hjust=0.5, size=14),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
          strip.background = element_rect(
            color="white", fill="white", linetype="solid"
          ),
          strip.text.x = element_text(size=12),
          strip.text.y = element_text(size=12)) +
    coord_fixed() + xlab(expression(paste("Cluster autocorrelation, ", r))) +
    ylab(expression(paste("Within-period ICC, ", rho))) +
    ggtitle(title, subtitle=stitle)
  
  rng <- ifelse(fullrange, "full", "restricted")
  corrstruct <- ifelse(corrtype==0, "BE", "DTD")
  rhouval <- ifelse((0 < rhou && rhou < 1), strsplit(as.character(rhou), "\\.")[[1]][2],
                    as.character(rhou))
  impper <- ifelse(imp, "_implementation", "")
  ggsave(paste0("plots/releff_SC_vs_SW_",
                corrstruct, "_", pereff, "_", rng, "_S_", S, "_rhou_",
                rhouval, impper, ".jpg"), p, width=6, height=5, units="in", dpi=800)
  ggsave(paste0("plots/releff_SC_vs_SW_",
                corrstruct, "_", pereff, "_", rng, "_S_", S, "_rhou_",
                rhouval, impper, ".pdf"), p, width=6, height=5, units="in", dpi=600)
  return(p)
}

releffSCSW_grid_multiplot_diffm_singleS <- function(S, K, corrtype, pereff,
                                                    fixedscale=FALSE, limits=c(0,2),
                                                    breaks=seq(0,2,0.2),
                                                    fullrange=FALSE, rhou=0, imp=FALSE){
  # Compare variances of complete SW designs and basic staircase designs with
  # larger cluster-period sizes (the same total number of participants), for a
  # range of correlation parameters, and for different cluster-period sizes
  # Inputs:
  #  S - number of unique treatment sequences
  #  K - number of clusters randomised to each sequence
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  #  fixedscale - Use particular limits and breaks for colour scale (T/F)
  #  limits - Range of values for colour scale shading
  #  breaks - Indices for colour scale legend and contour lines
  #  fullrange - Use full correlation parameter ranges (T/F)
  #  rhou - correlation between a participant's repeated measurements
  #          (default of 0 means each participant is measured just once)
  # Output:
  #  Contour plot of relative efficiencies (vartheta_SW/vartheta_SC)
  
  m1_SW <- 10
  m1_SC <- ((S+1)/2)*m1_SW
  m2_SW <- 100
  m2_SC <- ((S+1)/2)*m2_SW
  relvars_m1 <- gridvals_small(m1_SW, S, K, m1_SC, S, K,
                                corrtype, pereff, fullrange, rhou, imp)
  relvars_m1$m <- m1_SW
  relvars_m1$S <- S
  relvars_m2 <- gridvals_small(m2_SW, S, K, m2_SC, S, K,
                                corrtype, pereff, fullrange, rhou, imp)
  relvars_m2$m <- m2_SW
  relvars_m2$S <- S
  relvars <- bind_rows(
    relvars_m1,
    relvars_m2
  )
  
  if(fixedscale==TRUE){
    fillopt <- scale_fill_continuous_divergingx(name="Relative efficiency ",
                                                palette='RdBu', mid=1, rev=FALSE,
                                                limits=limits, breaks=breaks)
  }else{
    fillopt <- scale_fill_continuous_divergingx(name="Relative efficiency ",
                                                palette='RdBu', mid=1, rev=FALSE)
  }
  
  m.labs <- c(paste("mSW=10 \nmSC=", m1_SC), paste("mSW=100 \nmSC=", m2_SC))
  names(m.labs) <- c(m1_SW, m2_SW)
  S.labs <- paste("S = ", S)
  names(S.labs) <- S

  Kval <- ifelse(K==1, "", K)
  title <- bquote(paste("Relative efficiency, ",
                        var(hat(theta))[paste("SW(", .(S), ",", .(Kval), "K)")]/
                          var(hat(theta))[paste("SC(", .(S), ",", .(Kval), "K,1,1)")]))
  if(rhou==0 && imp){
    stitle <- "Including implementation periods"
  }else if(rhou!=0 && imp){
    stitle <- bquote(paste("Closed cohort with ", rho[u], "=", .(rhou),
                           ", including implementation periods"))
  }else if(rhou!=0 && isFALSE(imp)){
    stitle <- bquote(paste("Closed cohort with ", rho[u], "=", .(rhou)))
  }else{
    stitle <- NULL
  }
  
  p <- ggplot(relvars, aes(x=r, y=rho, z=releffSCSW)) + 
    geom_tile(aes(fill=releffSCSW)) +
    #    geom_contour_fill(breaks=breaks) +
    geom_contour2(color="black", breaks=breaks) +
    geom_text_contour(stroke=0.15, size=3, rotate=FALSE, breaks=breaks, skip=0,
                      label.placer=label_placer_fraction(0.5)) +
    facet_grid(
      m ~ S,
      labeller = labeller(m = m.labs, S=S.labs)
    ) +
    fillopt +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(aspect.ratio=3/8,
#          panel.spacing = unit(1.5, "lines"),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=14), legend.text=element_text(size=14),
          legend.position="bottom",
          plot.title=element_text(hjust=0.5, size=16),
          plot.subtitle=element_text(hjust=0.5, size=14),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
          strip.background = element_rect(
            color="white", fill="white", linetype="solid"
          ),
          strip.text.x = element_text(size=12),
          strip.text.y = element_text(size=12)) +
    coord_fixed() + xlab(expression(paste("Cluster autocorrelation, ", r))) +
    ylab(expression(paste("Within-period ICC, ", rho))) +
    ggtitle(title, subtitle=stitle)
  
  rng <- ifelse(fullrange, "full", "restricted")
  corrstruct <- ifelse(corrtype==0, "BE", "DTD")
  rhouval <- ifelse((0 < rhou && rhou < 1), strsplit(as.character(rhou), "\\.")[[1]][2],
                    as.character(rhou))
  impper <- ifelse(imp, "_implementation", "")
  ggsave(paste0("plots/releff_diffm_SC_vs_SW_",
                corrstruct, "_", pereff, "_", rng, "_S_", S, "_rhou_",
                rhouval, impper, ".jpg"), p, width=7, height=5, units="in", dpi=800)
  ggsave(paste0("plots/releff_diffm_SC_vs_SW_",
                corrstruct, "_", pereff, "_", rng, "_S_", S, "_rhou_",
                rhouval, impper, ".pdf"), p, width=7, height=5, units="in", dpi=600)
  return(p)
}

releff_SW_extendedSC_singleS <- function(S, K_SW, corrtype, pereff,
                                         fixedscale=FALSE, limits=c(0,2.5),
                                         breaks=seq(0,2.5,0.5),
                                         fullrange=FALSE, rhou=0, imp=FALSE){
  # Compare variances of complete SW designs and basic staircase designs with
  # more clusters (the same total number of participants), for a range of
  # correlation parameters, and for different cluster-period sizes
  # Inputs:
  #  S - number of unique treatment sequences
  #  K_SW - number of times each sequence is repeated for SW design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  #  fixedscale - Use particular limits and breaks for colour scale (T/F)
  #  limits - Range of values for colour scale shading
  #  breaks - Indices for colour scale legend and contour lines
  #  fullrange - Use full correlation parameter ranges (T/F)
  #  rhou - correlation between a participant's repeated measurements
  #          (default of 0 means each participant is measured just once)
  # Output:
  #  Contour plot of relative efficiencies (vartheta_SW/vartheta_SC)
  
  m1 <- 10
  m2 <- 100
  K_SC <- K_SW*(S+1)/2
  relvars_m1 <- gridvals_small(m1, S, K_SW, m1, S, K_SC,
                                corrtype, pereff, fullrange, rhou, imp)
  relvars_m1$m <- m1
  relvars_m1$S <- S
  relvars_m2 <- gridvals_small(m2, S, K_SW, m2, S, K_SC,
                                corrtype, pereff, fullrange, rhou, imp)
  relvars_m2$m <- m2
  relvars_m2$S <- S

  relvars <- bind_rows(
    relvars_m1,
    relvars_m2
  )
  
  if(fixedscale==TRUE){
    fillopt <- scale_fill_continuous_divergingx(name="Relative efficiency ",
                                                palette='RdBu', mid=1, rev=FALSE,
                                                limits=limits, breaks=breaks)
  }else{
    fillopt <- scale_fill_continuous_divergingx(name="Relative efficiency ",
                                                palette='RdBu', mid=1, rev=FALSE)
  }
  
  S.labs <- paste("S = ", S, ", q = ", K_SC)
  names(S.labs) <- S
  m.labs <- c(paste("m = ", m1), paste("m = ", m2))
  names(m.labs) <- c(m1, m2)
  
  KSW <- ifelse(K_SW==1, "", K_SW)
  title <- bquote(paste("Relative efficiency, ",
                        var(hat(theta))[paste("SW(", .(S), ",", .(KSW), "K)")]/
                          var(hat(theta))[paste("SC(", .(S), ",", .(K_SC), "K,1,1)")]))
  if(rhou==0 && imp){
    stitle <- "Including implementation periods"
  }else if(rhou!=0 && imp){
    stitle <- bquote(paste("Closed cohort with ", rho[u], "=", .(rhou),
                           ", including implementation periods"))
  }else if(rhou!=0 && isFALSE(imp)){
    stitle <- bquote(paste("Closed cohort with ", rho[u], "=", .(rhou)))
  }else{
    stitle <- NULL
  }
  
  p <- ggplot(relvars, aes(x=r, y=rho, z=releffSCSW)) +
    geom_tile(aes(fill=releffSCSW)) +
    #    geom_contour_fill(breaks=breaks) +
    geom_contour2(color="black", breaks=breaks) +
    geom_text_contour(stroke=0.15, size=3, rotate=FALSE, breaks=breaks, skip=0,
                      label.placer=label_placer_fraction(0.5)) +
    fillopt +
    facet_grid(
      m ~ S,
      labeller = labeller(m = m.labs, S = S.labs)
    ) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(aspect.ratio=3/8,
#          panel.spacing = unit(1.5, "lines"),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=14), legend.text=element_text(size=14),
          legend.position="bottom",
          plot.title=element_text(hjust=0.5, size=16),
          plot.subtitle=element_text(hjust=0.5, size=14),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
          strip.background = element_rect(
            color="white", fill="white", linetype="solid"
          ),
          strip.text.x = element_text(size=12),
          strip.text.y = element_text(size=12)) +
    coord_fixed() + xlab(expression(paste("Cluster autocorrelation, ", r))) +
    ylab(expression(paste("Within-period ICC, ", rho))) +
    ggtitle(title, subtitle=stitle)
  
  rng <- ifelse(fullrange, "full", "restricted")
  corrstruct <- ifelse(corrtype==0, "BE", "DTD")
  rhouval <- ifelse((0 < rhou && rhou < 1), strsplit(as.character(rhou), "\\.")[[1]][2],
                    as.character(rhou))
  impper <- ifelse(imp, "_implementation", "")
  ggsave(paste0("plots/releff_extendedSC_Sq11_vs_SW_S", K_SW, "_",
                corrstruct, "_", pereff, "_", rng, "_S_", S, "_rhou_",
                rhouval, impper, ".jpg"), p, width=7, height=5, units="in", dpi=800)
  ggsave(paste0("plots/releff_extendedSC_Sq11_vs_SW_S", K_SW, "_",
                corrstruct, "_", pereff, "_", rng, "_S_", S, "_rhou_",
                rhouval, impper, ".pdf"), p, width=7, height=5, units="in", dpi=600)
  return(p)
}

releffSC_SW_psi <- function(Svals, K_SW, m, rho0, pereff, diffK=FALSE){
  # Calculates the relative efficiency (variance of SW to SC), for all
  # combinations of psi and S values
  
  psivals <- seq(0, 0.995, 0.005)
  vals <- expand.grid(psi=psivals, S=Svals, K_SW=K_SW)

  if(diffK){
    vals <- vals %>%
      mutate(
        K_SC = (S+1)/2
      )
  }else{
    vals <- vals %>%
      mutate(
        K_SC = K_SW
      )
  }
    
  vals$varSW <- with(
    vals,
    sapply(1:nrow(vals), function(j){
      VarSW_alt_psi(psi[j], S[j], K_SW[j], m, rho0)
    }
    )
  )

  if(pereff=='cat'){
    vals$varSC <- with(
      vals,
      sapply(1:nrow(vals), function(j){
        VarSCcat_psi(psi[j], S[j], K_SC[j], m, rho0)
      }
      )
    )
  }else if(pereff=='lin'){
    vals$varSC <- with(
      vals,
      sapply(1:nrow(vals), function(j){
        VarSClin_psi(psi[j], S[j], K_SC[j], m, rho0)
      }
      )
    )
  }

  vals <- vals %>%
    mutate(
      releffSCvsSW = varSW/varSC
    )
  return(vals)
}

RE_lines_psi <- function(Svals, K_SW, pereff, diffK=FALSE){
  # Inputs:
  #  Svals - different numbers of sequences to plot together, e.g. c(3,6)
  #  K - number of clusters to randomise to each sequence
  #  m - number of subjects per cluster-period
  #  rho0 - within-period intracluster correlation
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  # Output:
  #  Plot of relative efficiencies for designs with different numbers of sequences
  # Assumptions:
  #  Total variance of 1
  #  Block-exchangeable intracluster correlation structure
  
  # Note: RE is independent of m and rho0 for comparisons of stepped wedge with
  #       embedded and extended staircase designs
  REs <- releffSC_SW_psi(Svals=Svals, K_SW=K_SW, m=10, rho0=0.05, pereff=pereff, diffK=diffK)
  
  Kval <- ifelse(K_SW==1, "", K_SW)
  if(length(Svals)==1 && diffK==FALSE){
    title <- bquote(paste("Relative efficiency, ",
                          var(hat(theta))[paste("SW(", .(Svals), ",", .(Kval), "K)")]/
                            var(hat(theta))[paste("SC(", .(Svals), ",", .(Kval), "K,1,1)")]))
  }else if(length(Svals)==1 && diffK==TRUE){
    KSC <- (S+1)/2
    title <- bquote(paste("Relative efficiency, ",
                          var(hat(theta))[paste("SW(", .(Svals), ",", .(Kval), "K)")]/
                            var(hat(theta))[paste("SC(", .(Svals), ",", .(KSC), "K,1,1)")]))
  }else if(length(Svals)!=1 && diffK==FALSE){
    title <- bquote(paste("Relative efficiency, ",
                          var(hat(theta))[paste("SW(S,", .(Kval), "K)")]/
                            var(hat(theta))[paste("SC(S,", .(Kval), "K,1,1)")]))
  }else{
    KSC <- "((S+1)/2)"
    title <- bquote(paste("Relative efficiency, ",
                          var(hat(theta))[paste("SW(S,", .(Kval), "K)")]/
                            var(hat(theta))[paste("SC(S,", .(KSC), "K,1,1)")]))
  }

  p <- ggplot(data=REs, aes(x=psi, y=releffSCvsSW, colour=as.factor(S), linetype=as.factor(S))) +
    geom_hline(yintercept=1, linewidth=1.5, color="gray75", linetype="solid") +
    geom_line(linewidth=2) +
    scale_color_manual(values=colorRampPalette(c("lightgreen", "darkgreen"))(length(Svals)),
                       name=expression(paste("Number of sequences, ", S))) +
    xlab(expression(paste("Correlation between cluster-period means, ", psi))) +
    ylab("Relative efficiency") +
    labs(title=title,
         colour=expression(paste("Number of sequences, ", S)),
         linetype=expression(paste("Number of sequences, ", S))) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=16),
          plot.subtitle=element_text(hjust=0.5, size=14),
          axis.title=element_text(size=12), axis.text=element_text(size=12),
          legend.key.width = unit(1.5, "cm"), legend.title=element_text(size=14),
          legend.text=element_text(size=14), legend.position="bottom")
  return(p)
}
