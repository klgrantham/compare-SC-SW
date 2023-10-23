# Generate paper plots

source('paper_plots_comparisons.R')

# Figure 2: Correlation between cluster-period means
p2 <- corr_grid_multiplot(fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))
ggsave(paste0("plots/figure2.eps"), p2, width=6, height=5, units="in", dpi=800)

corr_grid_multiplot(fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2), fullrange=TRUE)

# Figure 3: Weight to horizontal comparisons
p3 <- weight_grid_multiplot()
ggsave(paste0("plots/figure3.eps"), p3, width=7, height=4, units="in", dpi=800)

### Relative efficiency, limited parameter range ###

## Embedded staircase vs stepped wedge ##

# Figure 4: Block-exchangeable, categorical period effects
p4 <- releffSCSW_grid_multiplot_corr_diffS(3, 9, 1, 1, 1, 1, 0, 'cat', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))
ggsave(paste0("plots/figure4.eps"), p4, width=9, height=5, units="in", dpi=800)

# Figure S1: Discrete-time decay, categorical period effects
releffSCSW_grid_multiplot_corr_diffS(3, 9, 1, 1, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))

# Figure S4: Block-exchangeable, linear period effects
releffSCSW_grid_multiplot_corr_diffS(3, 9, 1, 1, 1, 1, 0, 'lin', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))

# Figure S7: Block-exchangeable, categorical period effects, full range
releffSCSW_grid_multiplot_corr_diffS(3, 9, 1, 1, 1, 1, 0, 'cat', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2), fullrange=TRUE)

## Staircase with larger cluster-period size vs stepped wedge ##

# Figure 5: Block-exchangeable, categorical period effects
p5 <- releffSCSW_grid_multiplot_diffm_diffS(3, 9, 1, 1, 0, 'cat', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.5))
ggsave(paste0("plots/figure5.eps"), p5, width=9, height=5, units="in", dpi=800)

# Figure S2: Discrete-time decay, categorical period effects
releffSCSW_grid_multiplot_diffm_diffS(3, 9, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.5))

# Figure S5: Block-exchangeable, linear period effects
releffSCSW_grid_multiplot_diffm_diffS(3, 9, 1, 1, 0, 'lin', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5))

# Figure S8: Block-exchangeable, categorical period effects, full range
releffSCSW_grid_multiplot_diffm_diffS(3, 9, 1, 1, 0, 'cat', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.5), fullrange=TRUE)

## Extended staircase vs stepped wedge ##

# Figure 6: Block-exchangeable, categorical period effects
p6 <- releff_SW_extendedSC_diffS(3, 9, 1, 1, 1, 0, 'cat', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5))
ggsave(paste0("plots/figure6.eps"), p6, width=9, height=5, units="in", dpi=800)

# Figure S3: Discrete-time decay, categorical period effects
releff_SW_extendedSC_diffS(3, 9, 1, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,5.0), breaks=seq(0,5.0,1.0))

# Figure S6: Block-exchangeable, linear period effects
releff_SW_extendedSC_diffS(3, 9, 1, 1, 1, 0, 'lin', fixedscale=TRUE, limits=c(0,3.0), breaks=seq(0,3.0,0.5))

# Figure S9: Block-exchangeable, categorical period effects, full range
releff_SW_extendedSC_diffS(3, 9, 1, 1, 1, 0, 'cat', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5), fullrange=TRUE)


## Trial examples

# PROMPT trial
SWvar <- CRTvartheta(20, SWdesmat(5, 8), 0.03, 1, 0, 'cat')
pow(SWvar, 0.15)
SCvar_m20 <- VarSCcat(20, 5, 8, 0.03, 1)
pow(SCvar_m20, 0.15)
releff_SW_SCm20 <- SWvar/SCvar_m20

SCvar_m30 <- VarSCcat(30, 5, 8, 0.03, 1)
pow(SCvar_m30, 0.15)
releff_SW_SCm30 <- SWvar/SCvar_m30

SWvar_dtd <- CRTvartheta(20, SWdesmat(5, 8), 0.032, 0.97, 1, 'cat')
SCvar_m20_dtd <- VarSCcat(20, 5, 8, 0.032, 0.97)
pow(SCvar_m20_dtd, 0.15)
releff_SW_SCm20_dtd <- SWvar_dtd/SCvar_m20_dtd

SCvar_m30_dtd <- VarSCcat(30, 5, 8, 0.032, 0.97)
pow(SCvar_m30_dtd, 0.15)
releff_SW_SCm30_dtd <- SWvar_dtd/SCvar_m30_dtd
