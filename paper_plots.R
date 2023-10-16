# Generate paper plots

source('paper_plots_comparisons.R')

# Correlation between cluster-period means
corr_grid_multiplot(fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))
corr_grid_multiplot(fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2), fullrange=TRUE)

# Weight to horizontal comparisons
weight_grid_multiplot()

# Scaled variances of stepped wedge and staircase designs for different
# parameter configurations
VarSCbasic_SW_multi_line_plot()

# Variances of stepped wedge and staircase designs for different parameter
# configurations
VarSCbasic_SW_multi_line_plot_a(3, 1)
VarSCbasic_SW_multi_line_plot_a(10, 1)

## Relative efficiency, limited parameter range

# Embedded staircase vs stepped wedge

# Block-exchangeable, categorical period effects
releffSCSW_grid_multiplot_corr_diffS(3, 9, 1, 1, 1, 1, 0, 'cat', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))

# Discrete-time decay, categorical period effects
releffSCSW_grid_multiplot_corr_diffS(3, 9, 1, 1, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))

# Block-exchangeable, linear period effects
releffSCSW_grid_multiplot_corr_diffS(3, 9, 1, 1, 1, 1, 0, 'lin', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))

# Discrete-time decay, linear period effects
releffSCSW_grid_multiplot_corr_diffS(3, 9, 1, 1, 1, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))

# Block-exchangeable, categorical period effects, full range
releffSCSW_grid_multiplot_corr_diffS(3, 9, 1, 1, 1, 1, 0, 'cat', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2), fullrange=TRUE)


#releffSCSW_grid_multiplot_corr(3, 1, 3, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2), genK=TRUE)
#releffSCSW_grid_multiplot_corr(3, 1, 3, 1, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2), genK=TRUE)
#releffSCSW_grid_multiplot_corr(10, 1, 10, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2), genK=TRUE)
#releffSCSW_grid_multiplot_corr(10, 1, 10, 1, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2), genK=TRUE)

# Staircase with larger cluster-period size vs stepped wedge

# Block-exchangeable, categorical period effects
releffSCSW_grid_multiplot_diffm_diffS(3, 9, 1, 1, 0, 'cat', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.5))

# Discrete-time decay, categorical period effects
releffSCSW_grid_multiplot_diffm_diffS(3, 9, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.5))

# Block-exchangeable, linear period effects
releffSCSW_grid_multiplot_diffm_diffS(3, 9, 1, 1, 0, 'lin', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5))

# Discrete-time decay, linear period effects
releffSCSW_grid_multiplot_diffm_diffS(3, 9, 1, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5))

# Block-exchangeable, categorical period effects, full range
releffSCSW_grid_multiplot_diffm_diffS(3, 9, 1, 1, 0, 'cat', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.5), fullrange=TRUE)


#releffSCSW_grid_multiplot_diffm(3, 1, 3, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5), genK=TRUE)
#releffSCSW_grid_multiplot_diffm(3, 1, 3, 1, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5), genK=TRUE)

#releffSCSW_grid_multiplot_diffm(10, 1, 10, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5), genK=TRUE)
#releffSCSW_grid_multiplot_diffm(10, 1, 10, 1, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5), genK=TRUE)

# Extended staircase vs stepped wedge

# Block-exchangeable, categorical period effects
releff_SW_extendedSC_diffS(3, 9, 1, 1, 1, 0, 'cat', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5))

# Discrete-time decay, categorical period effects
releff_SW_extendedSC_diffS(3, 9, 1, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,5.0), breaks=seq(0,5.0,1.0))

# Block-exchangeable, linear period effects
releff_SW_extendedSC_diffS(3, 9, 1, 1, 1, 0, 'lin', fixedscale=TRUE, limits=c(0,3.0), breaks=seq(0,3.0,0.5))

# Discrete-time decay, linear period effects
releff_SW_extendedSC_diffS(3, 9, 1, 1, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,5.0), breaks=seq(0,5.0,1.0))

# Block-exchangeable, categorical period effects
releff_SW_extendedSC_diffS(3, 9, 1, 1, 1, 0, 'cat', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5), fullrange=TRUE)


#releffSCSW_grid_multiplot_corr(3, 1, 3, 2, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.5))
#releffSCSW_grid_multiplot_corr(3, 1, 3, 2, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.5))


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
