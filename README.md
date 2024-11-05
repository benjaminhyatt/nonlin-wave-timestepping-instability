# nonlin-wave-timestepping-instability

This repository contains scripts used to generate the data and figures presented in an upcoming manuscript "Multiple scales analysis of a nonlinear timestepping instability in simulations of solitons" by Hyatt et al. Below is a description of the contents of the repository: 

## Soliton propagation initial value problems
`kdv_ivp.py` implements Dedalus version 3 to propagate solitons on a periodic interval using IMEX timestepping schemes to inspect their (in)stability. When the scheme is unstable we run the simulation all the way through blow-up, whereas if the scheme is stable (i.e., RK443) we specify a late stop time. Each simulation can be configured to yield a variety of outputs, such as snapshots of the solution and its (squared) L2 norm versus time. The cadence of outputs can be adjusted depending on the time scale of the simulation. 

## von Neumann analysis eigenvalue problems

## Processing scripts

## Processed data

## Plotting scripts
