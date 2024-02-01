# Generate paper plots

source('paper_plots_comparisons.R')

# Figure 2: Correlation between cluster-period means
p2 <- corr_grid_multiplot(fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))
ggsave(paste0("plots/figure2.jpg"), p2, width=6, height=5, units="in", dpi=600)
ggsave(paste0("plots/figure2.eps"), p2, width=6, height=5, units="in", dpi=800)

corr_grid_multiplot(fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2), fullrange=TRUE)

# Figure 3: Weight to horizontal comparisons
p3 <- weight_grid_multiplot()
ggsave(paste0("plots/figure3.jpg"), p3, width=7, height=4, units="in", dpi=600)
ggsave(paste0("plots/figure3.eps"), p3, width=7, height=4, units="in", dpi=800)

### Relative efficiency, limited parameter range ###

## Embedded staircase vs stepped wedge ##

# Figure 4: Block-exchangeable, categorical period effects
p4 <- releffSCSW_grid_multiplot_corr_diffS(S1=3, S2=9, K_SW=1, K_SC=1,
       corrtype=0, pereff='cat', fixedscale=TRUE, limits=c(0,1),
       breaks=seq(0,1,0.2), fullrange=FALSE, rhou=0, imp=FALSE)
ggsave(paste0("plots/figure4.jpg"), p4, width=9, height=5, units="in", dpi=600)
ggsave(paste0("plots/figure4.eps"), p4, width=9, height=5, units="in", dpi=800)

# Figure S1: Discrete-time decay, categorical period effects
pS1 <- releffSCSW_grid_multiplot_corr_diffS(S1=3, S2=9, K_SW=1, K_SC=1,
        corrtype=1, pereff='cat', fixedscale=TRUE, limits=c(0,1),
        breaks=seq(0,1,0.2), fullrange=FALSE, rhou=0, imp=FALSE)
ggsave(paste0("plots/figureS1.pdf"), pS1, width=9, height=5, units="in", dpi=600)

# Figure S4: Block-exchangeable, linear period effects
pS4 <- releffSCSW_grid_multiplot_corr_diffS(S1=3, S2=9, K_SW=1, K_SC=1,
        corrtype=0, pereff='lin', fixedscale=TRUE, limits=c(0,1),
        breaks=seq(0,1,0.2), fullrange=FALSE, rhou=0, imp=FALSE)
ggsave(paste0("plots/figureS4.pdf"), pS4, width=9, height=5, units="in", dpi=600)

# Figure S7: Block-exchangeable, categorical period effects, full range
pS7 <- releffSCSW_grid_multiplot_corr_diffS(S1=3, S2=9, K_SW=1, K_SC=1,
        corrtype=0, pereff='cat', fixedscale=TRUE, limits=c(0,1),
        breaks=seq(0,1,0.2), fullrange=TRUE, rhou=0, imp=FALSE)
ggsave(paste0("plots/figureS7.pdf"), pS7, width=9, height=5, units="in", dpi=600)

## Staircase with larger cluster-period size vs stepped wedge ##

# Figure 5: Block-exchangeable, categorical period effects
p5 <- releffSCSW_grid_multiplot_diffm_diffS(S1=3, S2=9, corrtype=0,
        pereff='cat', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.2),
        fullrange=FALSE, rhou=0, imp=FALSE)
ggsave(paste0("plots/figure5.jpg"), p5, width=9, height=5, units="in", dpi=600)
ggsave(paste0("plots/figure5.eps"), p5, width=9, height=5, units="in", dpi=800)

# Figure S2: Discrete-time decay, categorical period effects
pS2 <- releffSCSW_grid_multiplot_diffm_diffS(S1=3, S2=9, corrtype=1,
         pereff='cat', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.5),
         fullrange=FALSE, rhou=0, imp=FALSE)
ggsave(paste0("plots/figureS2.pdf"), pS2, width=9, height=5, units="in", dpi=600)

# Figure S5: Block-exchangeable, linear period effects
pS5 <- releffSCSW_grid_multiplot_diffm_diffS(S1=3, S2=9, corrtype=0,
         pereff='lin', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5),
         fullrange=FALSE, rhou=0, imp=FALSE)
ggsave(paste0("plots/figureS5.pdf"), pS5, width=9, height=5, units="in", dpi=600)

# Figure S8: Block-exchangeable, categorical period effects, full range
pS8 <- releffSCSW_grid_multiplot_diffm_diffS(S1=3, S2=9, corrtype=0,
         pereff='cat', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.5),
         fullrange=TRUE, rhou=0, imp=FALSE)
ggsave(paste0("plots/figureS8.pdf"), pS8, width=9, height=5, units="in", dpi=600)

## Extended staircase vs stepped wedge ##

# Figure 6: Block-exchangeable, categorical period effects
p6 <- releff_SW_extendedSC_diffS(S1=3, S2=9, K_SW=1, corrtype=0, pereff='cat',
       fixedscale=TRUE, limits=c(0,2.4), breaks=seq(0,2.4,0.2), fullrange=FALSE,
       rhou=0, imp=FALSE)
ggsave(paste0("plots/figure6.jpg"), p6, width=9, height=5, units="in", dpi=600)
ggsave(paste0("plots/figure6.eps"), p6, width=9, height=5, units="in", dpi=800)

# Figure S3: Discrete-time decay, categorical period effects
pS3 <- releff_SW_extendedSC_diffS(S1=3, S2=9, K_SW=1, corrtype=1, pereff='cat',
        fixedscale=TRUE, limits=c(0,5.0), breaks=seq(0,5.0,1.0), fullrange=FALSE,
        rhou=0, imp=FALSE)
ggsave(paste0("plots/figureS3.pdf"), pS3, width=9, height=5, units="in", dpi=600)

# Figure S6: Block-exchangeable, linear period effects
pS6 <- releff_SW_extendedSC_diffS(S1=3, S2=9, K_SW=1, corrtype=0, pereff='lin',
        fixedscale=TRUE, limits=c(0,3.0), breaks=seq(0,3.0,0.5), fullrange=FALSE,
        rhou=0, imp=FALSE)
ggsave(paste0("plots/figureS6.pdf"), pS6, width=9, height=5, units="in", dpi=600)

# Figure S9: Block-exchangeable, categorical period effects, full range
pS9 <- releff_SW_extendedSC_diffS(S1=3, S2=9, K_SW=1, corrtype=0, pereff='cat',
        fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5), fullrange=TRUE,
        rhou=0, imp=FALSE)
ggsave(paste0("plots/figureS9.pdf"), pS9, width=9, height=5, units="in", dpi=600)


## Trial examples

# PROMPT trial
SWvar <- CRTvartheta(m=20, Xmat=SWdesmat(5, 8), rho0=0.03, r=1, corrtype=0,
                     pereff='cat', rhou=0)
pow(SWvar, 0.15)
SCvar_m20 <- VarSCcat(m=20, S=5, K=8, rho0=0.03, r=1, rhou=0)
pow(SCvar_m20, 0.15)
releff_SW_SCm20 <- SWvar/SCvar_m20

SCvar_m30 <- VarSCcat(m=30, S=5, K=8, rho0=0.03, r=1, rhou=0)
pow(SCvar_m30, 0.15)
releff_SW_SCm30 <- SWvar/SCvar_m30

SWvar_dtd <- CRTvartheta(m=20, Xmat=SWdesmat(5, 8), rho0=0.032, r=0.97,
                         corrtype=1, pereff='cat', rhou=0)
SCvar_m20_dtd <- VarSCcat(m=20, S=5, K=8, rho0=0.032, r=0.97, rhou=0)
pow(SCvar_m20_dtd, 0.15)
releff_SW_SCm20_dtd <- SWvar_dtd/SCvar_m20_dtd

SCvar_m30_dtd <- VarSCcat(m=30, S=5, K=8, rho0=0.032, r=0.97, rhou=0)
pow(SCvar_m30_dtd, 0.15)
releff_SW_SCm30_dtd <- SWvar_dtd/SCvar_m30_dtd
