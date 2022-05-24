using DifferentialEquations
using Plots
using LinearAlgebra
using Distributions
using Random

Random.seed!(123)
# Initial condition
d0 = MvNormal(ones(2), Diagonal(ones(2)))
x0 = rand(d0)

# System
f(x) = [-x[1] + 2*x[2];
        -0.05*x[2]^3 - 3*x[1] - x[2]]
B = [0.5, -2.0]

# (Unknown) Diffusion Parameter 
sigma(t, x) = 0.01*x[2]^2
# Sim time
tspan = (0., 5.)

# Dual of the control contraction metric matrix (or as a function if state-dependent)
W = [4.26662 -0.921495; -0.921495 3.73143]
λ = 1.74

# Define desired regulation point or trajectories (as functions)
xs = [0.0, 0.0]
us = 0.0

# Uncertainty and upper bound
h(t, x) = 2*sin(2*t)
Δh = 2.0
Sigma(t, x) = 0.01*x[1]^2
Δg = 0.04

# L1 filter bandwidth and adaptation rate
ω = 90 #500 for another case
Γ = 4e7 #40 for another case

sys_p = sys_params(f, B)
ccm_p = ccm_params(xs, us, λ, W)
l1_p = l1_params(ω, Γ, Δh, Δg)

# L1 Reference sysem with perturbations AKA (non-implementable) ideal
# performance of the L1 system with full knowledge of the uncertainty
ref_sys = reference_system(sys_p, ccm_p, ω, h, sigma, Sigma)
ref_sol = solve(ref_sys, xs0, tspan, LambaEM(), progress = true, progress_steps = 1, saveat = 0.1, maxiters=5e7)

# L1 + CCM system to handle disturbances
l1_sys = l1_system(sys_p, ccm_p, l1_p, h, sigma, Sigma)
l1_sol = solve(l1_sys, x0, tspan, LambaEM(), progress = true, progress_steps = 1, saveat = 0.1, maxiters=5e7)

# Generate trajectory points
ref_sol_x1=Vector{Float64}();
ref_sol_x2=Vector{Float64}();
l1_sol_x1=Vector{Float64}();
l1_sol_x2=Vector{Float64}();
for i=1:51
        append!(l1_sol_x1,l1_sol[i][1])
        append!(l1_sol_x2,l1_sol[i][2])
        append!(ref_sol_x1,ref_sol[i][1])
        append!(ref_sol_x2,ref_sol[i][2])
end

#MSE between actual system and reference system
mseactualandref=(l1_sol_x1-ref_sol_x1).^2+(l1_sol_x2-ref_sol_x2).^2;

#MSE between reference system and ideal system
mserefandideal=(ref_sol_x1).^2+(ref_sol_x2).^2;

#MSE between actual system and ideal system
mseactualandideal=(l1_sol_x1).^2+(l1_sol_x2).^2;

#2-Wasserstein Distance between actual and ideal system
Wass2actualandideal=mseactualandideal.^0.5;

#Generate Plot
t=0:0.1:5;
p1=plot(t,mseactualandideal)
p2=plot(t,mseactualandref)
p3=plot(t,mserefandideal)
p4=plot(t,Wass2actualandideal)
