# nonlin-wave-timestepping-instability

This repository contains scripts used to generate the data and figures presented in an upcoming manuscript "Multiple scales analysis of a nonlinear timestepping instability in simulations of solitons" by Hyatt et al. Below is a description of the contents of the repository: 

# Folders:
## `processed_data`
This folder consists of data that has been processed (using the scripts inside of `processing`) and is ready to be plotted. We provide the processed data rather than raw data due to the large volumes of raw data used in this work. 

## `plotting`
This folder consists of scripts that read in data from `processed_data` and output figures in `.eps` format. Scripts should be run from the main directory such that the `processed_data` folder is visible. 

## `processing`
This folder consists of scripts used to process data produced by any of the initial value problems, eigenvalue problems, and/or multiple scales analysis predictions described in this work. Due to the large volume of data produced in this work, we provide the scripts to generate your own data (see descriptions below). 

# Scripts in the main directory:  
## Soliton propagation initial value problems
`kdv_ivp.py` implements the IVP class in Dedalus version 3 (see [https://github.com/DedalusProject/dedalus](https://github.com/DedalusProject/dedalus)) to propagate solitons on a periodic interval using IMEX timestepping schemes to inspect their (in)stability. When the scheme is unstable we run the simulation all the way through blow-up, whereas if the scheme is stable (i.e., RK443) we specify a late stop time. Each simulation can be configured to yield a variety of outputs, such as snapshots of the solution and its (squared) L2 norm versus time. The cadence of outputs can be adjusted depending on the time scale of the simulation. 

## von Neumann analysis eigenvalue problems
`kdv_vn_evp.py` implements the EVP class in Dedalus version 3 to investigate whether an IMEX scheme will stably timestep perturbations to a background soliton under the linearized KdV equation. This script consists of functions for each of the six IMEX schemes considered in this work, and the setup can be generalized to other multi-step and multi-stage schemes. We provide functions implementing both `dense` and `sparse` strategies for solving the EVP which, for instance, may be used to investigate the full spectrum or just obtain information about the fastest growing modes, respectively.
NOTE: we also provide `solvers.py` which modifies the script of the same name found in the Dedalus version 3 repository [https://github.com/DedalusProject/dedalus/blob/master/dedalus/core/solvers.py](https://github.com/DedalusProject/dedalus/blob/master/dedalus/core/solvers.py). This is done in order to apply a phase shift operator after construction of the (unshifted) problem matrices, which is done in a slightly different manner depending on whether the scheme is multi-step or multi-stage.

## Multiple scales analysis: numerical integration of the nontrivial solvability conditions
`ms_integrate_predictions.py` initiates the solution of the nontrivial solvability conditions based on those derived using the multiple scales framework developed in this work. This script relies on `solvability_expr_final.py` (which stores the nonlinear ODEs for the effective soliton parameter c) and `kdv_integrate_finite_final.py` which uses the `scipy.integrate` library to numerically integrate the ODEs. 
