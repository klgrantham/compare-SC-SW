# Generate paper plots

source('paper_plots_comparisons.R')

# Figure 2: Correlation between cluster-period means
p2 <- corr_grid_multiplot(fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))
ggsave(paste0("plots/figure2.jpg"), p2, width=6, height=5, units="in", dpi=600)
ggsave(paste0("plots/figure2.eps"), p2, width=6, height=5, units="in", dpi=800)

corr_grid_multiplot(fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2), fullrange=TRUE)

### Relative efficiency ###

## Embedded staircase vs stepped wedge ##

# Figure 3: RE vs psi, block-exchangeable, categorical period effects
p3 <- RE_lines_psi(Svals=c(3,5,7,9,11), K_SW=1, pereff='cat', diffK=FALSE)
ggsave(paste0("plots/figure3.jpg"), p3, width=9, height=5, units="in", dpi=600)
ggsave(paste0("plots/figure3.eps"), p3, width=9, height=5, units="in", dpi=800)

# Figure 4: Contour plot, block-exchangeable, categorical period effects
p4 <- releffSCSW_grid_multiplot_corr_diffS(S1=3, S2=9, K_SW=1, K_SC=1,
       corrtype=0, pereff='cat', fixedscale=TRUE, limits=c(0,1),
       breaks=seq(0,1,0.2), fullrange=FALSE, rhou=0, imp=FALSE)
ggsave(paste0("plots/figure4.jpg"), p4, width=9, height=5, units="in", dpi=600)
ggsave(paste0("plots/figure4.eps"), p4, width=9, height=5, units="in", dpi=800)

# Figure S5: Contour plot, discrete-time decay, categorical period effects
pS5 <- releffSCSW_grid_multiplot_corr_diffS(S1=3, S2=9, K_SW=1, K_SC=1,
        corrtype=1, pereff='cat', fixedscale=TRUE, limits=c(0,1),
        breaks=seq(0,1,0.2), fullrange=FALSE, rhou=0, imp=FALSE)
ggsave(paste0("plots/figureS5.pdf"), pS5, width=9, height=5, units="in", dpi=600)

# Figure S8: Contour plot, block-exchangeable, linear period effects
pS8 <- releffSCSW_grid_multiplot_corr_diffS(S1=3, S2=9, K_SW=1, K_SC=1,
        corrtype=0, pereff='lin', fixedscale=TRUE, limits=c(0,1),
        breaks=seq(0,1,0.2), fullrange=FALSE, rhou=0, imp=FALSE)
ggsave(paste0("plots/figureS8.pdf"), pS8, width=9, height=5, units="in", dpi=600)

## Extended staircase vs stepped wedge ##

# Figure 5: RE vs psi, block-exchangeable, categorical period effects, same total number of participants
p5 <- RE_lines_psi(Svals=c(3,5,7,9,11), K_SW=1, pereff='cat', diffK=TRUE, same_tot=TRUE)
ggsave(paste0("plots/figure5.jpg"), p5, width=9, height=5, units="in", dpi=600)
ggsave(paste0("plots/figure5.eps"), p5, width=9, height=5, units="in", dpi=800)

# Figure 6: Contour plot, block-exchangeable, categorical period effects
p6 <- releff_SW_extendedSC_diffS(S1=3, S2=9, K_SW=1, corrtype=0, pereff='cat',
                                 fixedscale=TRUE, limits=c(0,2.4), breaks=seq(0,2.4,0.2), fullrange=FALSE,
                                 rhou=0, imp=FALSE)
ggsave(paste0("plots/figure6.jpg"), p6, width=9, height=5, units="in", dpi=600)
ggsave(paste0("plots/figure6.eps"), p6, width=9, height=5, units="in", dpi=800)

# Figure S1: RE vs psi, block-exchangeable, categorical period effects, fewer total participants
pS1 <- RE_lines_psi(Svals=c(3,5,7,9,11), K_SW=1, pereff='cat', diffK=TRUE, same_tot=FALSE)
ggsave(paste0("plots/figureS1.jpg"), pS1, width=9, height=5, units="in", dpi=600)
ggsave(paste0("plots/figureS1.pdf"), pS1, width=9, height=5, units="in", dpi=600)

# Figure S2: Contour plot, block-exchangeable, categorical period effects, fewer total participants
pS2 <- releff_SW_extendedSC_diffS(S1=5, S2=9, K_SW=1, corrtype=0, pereff='cat',
           fixedscale=TRUE, limits=c(0,1.4), breaks=seq(0,1.4,0.2), fullrange=FALSE,
           rhou=0, imp=FALSE, same_tot=FALSE)
ggsave(paste0("plots/figureS2.jpg"), pS2, width=9, height=5, units="in", dpi=600)
ggsave(paste0("plots/figureS2.pdf"), pS2, width=9, height=5, units="in", dpi=600)

# Figure S6: Contour plot, discrete-time decay, categorical period effects
pS6 <- releff_SW_extendedSC_diffS(S1=3, S2=9, K_SW=1, corrtype=1, pereff='cat',
                                  fixedscale=TRUE, limits=c(0,5.0), breaks=seq(0,5.0,1.0), fullrange=FALSE,
                                  rhou=0, imp=FALSE)
ggsave(paste0("plots/figureS6.pdf"), pS6, width=9, height=5, units="in", dpi=600)

# Figure S9: Contour plot, block-exchangeable, linear period effects
pS9 <- releff_SW_extendedSC_diffS(S1=3, S2=9, K_SW=1, corrtype=0, pereff='lin',
                                  fixedscale=TRUE, limits=c(0,3.0), breaks=seq(0,3.0,0.5), fullrange=FALSE,
                                  rhou=0, imp=FALSE)
ggsave(paste0("plots/figureS9.pdf"), pS9, width=9, height=5, units="in", dpi=600)


## Staircase with larger cluster-period size vs stepped wedge ##

# Figure 7: Contour plot, block-exchangeable, categorical period effects, same total number of participants
p7 <- releffSCSW_grid_multiplot_diffm_diffS(S1=3, S2=9, corrtype=0,
        pereff='cat', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.2),
        fullrange=FALSE, rhou=0, imp=FALSE)
ggsave(paste0("plots/figure7.jpg"), p7, width=9, height=5, units="in", dpi=600)
ggsave(paste0("plots/figure7.eps"), p7, width=9, height=5, units="in", dpi=800)

# Figure S3: Contour plot, block-exchangeable, categorical period effects, 25% increase
pS3 <- releffSCSW_grid_multiplot_diffm_diffS(S1=3, S2=9, corrtype=0,
          pereff='cat', fixedscale=TRUE, limits=c(0,1.0), breaks=seq(0,1.0,0.2),
          fullrange=FALSE, rhou=0, imp=FALSE, fixed_incr=TRUE, factor_incr=1.25)
ggsave(paste0("plots/figureS3.jpg"), pS3, width=9, height=5, units="in", dpi=600)
ggsave(paste0("plots/figureS3.pdf"), pS3, width=9, height=5, units="in", dpi=600)

# Figure S4: Contour plot, block-exchangeable, categorical period effects, 50% increase
pS4 <- releffSCSW_grid_multiplot_diffm_diffS(S1=3, S2=9, corrtype=0,
         pereff='cat', fixedscale=TRUE, limits=c(0,1.2), breaks=seq(0,1.2,0.2),
         fullrange=FALSE, rhou=0, imp=FALSE, fixed_incr=TRUE, factor_incr=1.5)
ggsave(paste0("plots/figureS4.jpg"), pS4, width=9, height=5, units="in", dpi=600)
ggsave(paste0("plots/figureS4.pdf"), pS4, width=9, height=5, units="in", dpi=600)

# Figure S7: Contour plot, discrete-time decay, categorical period effects
pS7 <- releffSCSW_grid_multiplot_diffm_diffS(S1=3, S2=9, corrtype=1,
         pereff='cat', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.5),
         fullrange=FALSE, rhou=0, imp=FALSE)
ggsave(paste0("plots/figureS7.pdf"), pS7, width=9, height=5, units="in", dpi=600)

# Figure S10: Contour plot, block-exchangeable, linear period effects
pS10 <- releffSCSW_grid_multiplot_diffm_diffS(S1=3, S2=9, corrtype=0,
         pereff='lin', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5),
         fullrange=FALSE, rhou=0, imp=FALSE)
ggsave(paste0("plots/figureS10.pdf"), pS10, width=9, height=5, units="in", dpi=600)


## Staircase vs stepped wedge, cohort designs ##

# Figure S11: Contour plot, block-exchangeable, categorical period effects, cohort w/ rhou=0.3
pS11 <- releffSCSW_grid_multiplot_corr_diffS(S1=3, S2=9, K_SW=1, K_SC=1,
         corrtype=0, pereff='cat', fixedscale=TRUE, limits=c(0,1),
         breaks=seq(0,1,0.2), fullrange=FALSE, rhou=0.3, imp=FALSE)
ggsave(paste0("plots/figureS11.pdf"), pS11, width=9, height=5, units="in", dpi=600)

# Figure S12: Contour plot, block-exchangeable, categorical period effects, cohort w/ rhou=0.7
pS12 <- releffSCSW_grid_multiplot_corr_diffS(S1=3, S2=9, K_SW=1, K_SC=1,
         corrtype=0, pereff='cat', fixedscale=TRUE, limits=c(0,1),
         breaks=seq(0,1,0.2), fullrange=FALSE, rhou=0.7, imp=FALSE)
ggsave(paste0("plots/figureS12.pdf"), pS12, width=9, height=5, units="in", dpi=600)

# Figure S13: Contour plot, block-exchangeable, categorical period effects, cohort w/ rhou=0.3
pS13 <- releff_SW_extendedSC_diffS(S1=3, S2=9, K_SW=1, corrtype=0, pereff='cat',
       fixedscale=TRUE, limits=c(0,2.4), breaks=seq(0,2.4,0.2), fullrange=FALSE,
       rhou=0.3, imp=FALSE)
ggsave(paste0("plots/figureS13.pdf"), pS13, width=9, height=5, units="in", dpi=600)

# Figure S14: Contour plot, block-exchangeable, categorical period effects, cohort w/ rhou=0.7
pS14 <- releff_SW_extendedSC_diffS(S1=3, S2=9, K_SW=1, corrtype=0, pereff='cat',
          fixedscale=TRUE, limits=c(0,2.4), breaks=seq(0,2.4,0.2), fullrange=FALSE,
          rhou=0.7, imp=FALSE)
ggsave(paste0("plots/figureS14.pdf"), pS14, width=9, height=5, units="in", dpi=600)

# Figure S15: Contour plot, block-exchangeable, categorical period effects, cohort w/ rhou=0.3
pS15 <- releffSCSW_grid_multiplot_diffm_diffS(S1=3, S2=9, corrtype=0,
                                             pereff='cat', fixedscale=TRUE, limits=c(0,2.2), breaks=seq(0,2.2,0.2),
                                             fullrange=FALSE, rhou=0.3, imp=FALSE)
ggsave(paste0("plots/figureS15.pdf"), pS15, width=9, height=5, units="in", dpi=600)

# Figure S16: Contour plot, block-exchangeable, categorical period effects, cohort w/ rhou=0.7
pS16 <- releffSCSW_grid_multiplot_diffm_diffS(S1=3, S2=9, corrtype=0,
                                              pereff='cat', fixedscale=TRUE, limits=c(0,2.2), breaks=seq(0,2.2,0.2),
                                              fullrange=FALSE, rhou=0.7, imp=FALSE)
ggsave(paste0("plots/figureS16.pdf"), pS16, width=9, height=5, units="in", dpi=600)


## Staircase vs stepped wedge, with implementation periods ##

# Figure S18: Embedded staircase vs. stepped wedge, block-exchangeable,
#  categorical period effects, implementation periods
# Note: We only consider 9-sequence designs here since we would not recommend a
#       3-sequence design with implementation periods
pS18 <- releffSCSW_grid_multiplot_corr_singleS(S=9, K=1, corrtype=0, pereff='cat',
         fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2), fullrange=FALSE,
         rhou=0, imp=TRUE)
ggsave(paste0("plots/figureS18.pdf"), pS18, width=6, height=5, units="in", dpi=600)

# Figure S19: Staircase with more clusters vs. stepped wedge,
#  block-exchangeable, categorical period effects, implementation periods
pS19 <- releff_SW_extendedSC_singleS(S=9, K_SW=1, corrtype=0, pereff='cat',
         fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5),
         fullrange=FALSE, rhou=0, imp=TRUE)
ggsave(paste0("plots/figureS19.pdf"), pS19, width=6, height=5, units="in", dpi=600)

# Figure S20: Staircase with larger cluster-period size vs. stepped wedge,
#  block-exchangeable, categorical period effects, implementation periods
pS20 <- releffSCSW_grid_multiplot_diffm_singleS(S=9, K=1, corrtype=0,
          pereff='cat', fixedscale=TRUE, limits=c(0,2), breaks=seq(0,2,0.5),
          fullrange=FALSE, rhou=0, imp=TRUE)
ggsave(paste0("plots/figureS20.pdf"), pS20, width=6, height=5, units="in", dpi=600)


## Trial example

# PROMPT trial
SWvar <- CRTvartheta(m=20, Xmat=SWdesmat(5, 8), rho0=0.032, r=0.93, corrtype=0,
                     pereff='cat', rhou=0)
pow(SWvar, 0.15)
SCvar_m20 <- VarSCcat(m=20, S=5, K=8, rho0=0.032, r=0.93, rhou=0)
pow(SCvar_m20, 0.15)
releff_SW_SCm20 <- SWvar/SCvar_m20

SCvar_m30 <- VarSCcat(m=30, S=5, K=8, rho0=0.032, r=0.93, rhou=0)
pow(SCvar_m30, 0.15)
releff_SW_SCm30 <- SWvar/SCvar_m30

SCvar_m20_50clust <- VarSCcat(m=20, S=5, K=10, rho0=0.032, r=0.93, rhou=0)
pow(SCvar_m20_50clust, 0.15)
releff_SW_SCm20_50clust <- SWvar/SCvar_m20_50clust

SCvar_m20_60clust <- VarSCcat(m=20, S=5, K=12, rho0=0.032, r=0.93, rhou=0)
pow(SCvar_m20_60clust, 0.15)
releff_SW_SCm20_60clust <- SWvar/SCvar_m20_60clust
