# StochasticMotionPlanning.jl

_This package is under active development_

Design safe controllers for nonlinear systems in the presence of disturbances and noises. The controllers ensure that the system state remains in a tunable safe bound around the desired state. Find out more about it in our paper:


Run all files in "src" folder before running files in "examples" folder, the correct order is:
1.  chebyshev.jl
2.  types.jl
3.  sys.jl
4.  ccm.jl
5.  l1.jl

Before running ex2.jl, run all files in "iterative LQG" folder. Note that all files in "iterative LQG" folder are running in Matlab.
