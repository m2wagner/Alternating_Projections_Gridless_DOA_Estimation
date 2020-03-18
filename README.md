# Alternating_Projections_Gridless_DOA_Estimation
MATLAB code for Alternating Projections based gridless direction of arrival (DOA) estimation

run_APG_arb.m: 
The main script from which every DOA estimation algorithm can be tested

ADMM_gridless.m:
Alternating direction method of multipliers (ADMM) for gridless DOA estimation (ULA only)

ADMM_NUA.m:
ADMM adapted for non-uniform arrays

AP_Gridless.m:
Gridless DOA estimation using Alternating Projections 

CBF_spec.m:
Conventional beamforming (CBF) spectrum

do_LASSO.m:
Least absolute shrinkage and selector operator (LASSO) for gridded DOA estimation

do_SBLML3.m:
Sparse Bayesian Learning (SBL) for gridded DOA estimation

find_CBF_peaks.m:
Finds the peaks of the CBF spectrum

gen_DOAs.m:
Randomly generated DOAs with specified wrap around distance separatation

gen_signals_SNR.m:
Generates synthetic array measurements according to specified signal to noise ratio (SNR)

irregular_rootMUSIC.m:
Irregular root-MUSIC algorithm for non-uniform arrays

ITP.m:
Irregular Toeplitz projection, projects a matrix to nearest specified irregular Toeplitz matrix

IVD.m:
Performs irregular Vandermonde decomposition

MS_roots.m:
Finds roots of null spectrum using manifold separation technique

null_spec_polynomial.m:
Outputs value of null spectrum polynomial at specified location

PSD.m:
Projection onto positive semi-definite cone (PSD)

rootMUSIC.m:
Performs root-MUSIC algorithm

Tproj.m:
Projection to nearest Toeplitz matrix

plot_NSF_NUA.m:
Plots null spectrum function for example non-uniform array

plot_NSF_ULA.m:
Plots null spectrum function for example uniform array

