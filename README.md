# nonlin-wave-timestepping-instability

This repository contains scripts used to generate the data and figures presented in an upcoming manuscript "Multiple scales analysis of a nonlinear timestepping instability in simulations of solitons" by Hyatt et al. Below is a description of the contents of the repository: 

## Soliton propagation initial value problems
`kdv_ivp.py` implements the IVP class in Dedalus version 3 (see [https://github.com/DedalusProject/dedalus](https://github.com/DedalusProject/dedalus)) to propagate solitons on a periodic interval using IMEX timestepping schemes to inspect their (in)stability. When the scheme is unstable we run the simulation all the way through blow-up, whereas if the scheme is stable (i.e., RK443) we specify a late stop time. Each simulation can be configured to yield a variety of outputs, such as snapshots of the solution and its (squared) L2 norm versus time. The cadence of outputs can be adjusted depending on the time scale of the simulation. 

## von Neumann analysis eigenvalue problems
`kdv_vn_evp.py` implements the EVP class in Dedalus version 3 to investigate whether an IMEX scheme will stably timestep perturbations to a background soliton under the linearized KdV equation. This script consists of functions for each of the six IMEX schemes considered in this work, and the setup can be generalized to other multi-step and multi-stage schemes. We provide functions implementing both `dense` and `sparse` strategies for solving the EVP which, for instance, may be used to investigate the full spectrum or just obtain information about the fastest growing modes, respectively.
NOTE: we also provide `solvers.py` which modifies the script of the same name found in the Dedalus version 3 repository [https://github.com/DedalusProject/dedalus/blob/master/dedalus/core/solvers.py](https://github.com/DedalusProject/dedalus/blob/master/dedalus/core/solvers.py). This is done in order to apply a phase shift operator after construction of the (unshifted) problem matrices, which is done in a slightly different manner depending on whether the scheme is multi-step or multi-stage.

## Processing scripts

## Processed data

## Plotting scripts
