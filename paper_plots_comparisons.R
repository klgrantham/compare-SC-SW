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
    fillopt <- scale_fill_distiller(name="Correlation", 
                                    palette="YlOrRd", direction=1,
                                    limits=limits, breaks=breaks)
  }else{
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
