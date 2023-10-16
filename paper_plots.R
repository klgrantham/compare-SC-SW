# Generate paper plots

source('paper_plots_comparisons.R')

# Figure 2: Correlation between cluster-period means
corr_grid_multiplot(fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))
corr_grid_multiplot(fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2), fullrange=TRUE)

# Figure 3: Weight to horizontal comparisons
weight_grid_multiplot()

### Relative efficiency, limited parameter range ###

## Embedded staircase vs stepped wedge ##

# Figure 4: Block-exchangeable, categorical period effects
releffSCSW_grid_multiplot_corr_diffS(3, 9, 1, 1, 1, 1, 0, 'cat', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))

# Figure S1: Discrete-time decay, categorical period effects
releffSCSW_grid_multiplot_corr_diffS(3, 9, 1, 1, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))

# Figure S4: Block-exchangeable, linear period effects
releffSCSW_grid_multiplot_corr_diffS(3, 9, 1, 1, 1, 1, 0, 'lin', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))

# Figure S7: Block-exchangeable, categorical period effects, full range
releffSCSW_grid_multiplot_corr_diffS(3, 9, 1, 1, 1, 1, 0, 'cat', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2), fullrange=TRUE)

## Staircase with larger cluster-period size vs stepped wedge ##

# Figure 5: Block-exchangeable, categorical period effects
releffSCSW_grid_multiplot_diffm_diffS(3, 9, 1, 1, 0, 'cat', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.5))

# Figure S2: Discrete-time decay, categorical period effects
releffSCSW_grid_multiplot_diffm_diffS(3, 9, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.5))

# Figure S5: Block-exchangeable, linear period effects
releffSCSW_grid_multiplot_diffm_diffS(3, 9, 1, 1, 0, 'lin', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5))

# Figure S8: Block-exchangeable, categorical period effects, full range
releffSCSW_grid_multiplot_diffm_diffS(3, 9, 1, 1, 0, 'cat', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.5), fullrange=TRUE)

## Extended staircase vs stepped wedge ##

# Figure 6: Block-exchangeable, categorical period effects
releff_SW_extendedSC_diffS(3, 9, 1, 1, 1, 0, 'cat', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5))

# Figure S3: Discrete-time decay, categorical period effects
releff_SW_extendedSC_diffS(3, 9, 1, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,5.0), breaks=seq(0,5.0,1.0))

# Figure S6: Block-exchangeable, linear period effects
releff_SW_extendedSC_diffS(3, 9, 1, 1, 1, 0, 'lin', fixedscale=TRUE, limits=c(0,3.0), breaks=seq(0,3.0,0.5))

# Figure S9: Block-exchangeable, categorical period effects, full range
releff_SW_extendedSC_diffS(3, 9, 1, 1, 1, 0, 'cat', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5), fullrange=TRUE)


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
