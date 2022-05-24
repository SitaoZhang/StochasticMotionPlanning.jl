using DifferentialEquations
using Plots
using LinearAlgebra
using MAT
using Random, Distributions
using LaTeXStrings

# Initial condition
x0 = [-10; 10]
xs0=[-10; 10]

# System
f(x) = [-x[1] + 2*x[2];
        -0.05*x[2]^3 - 3*x[1] - x[2]]
B = [0.5, -2.0]

# (Unknown) Diffusion Parameter 
Sigma(t, x) = 0.01*(x[2])^2;


# Dual of the control contraction metric matrix (or as a function if state-dependent)
W = [4.26662 -0.921495; -0.921495 3.73143]
λ = 1.74


# Uncertainty and upper bound
h(t, x) = 2*sin(2*t)
Δh = 2.0
sigma(t, x) = 0.01*x[1]^2
Δg = 2.0

# L1 filter bandwidth and adaptation rate
ω = 50000 #90 for another case
Γ = 4e9 #4e4 for another case

sys_p = sys_params(f, B)

l1_p = l1_params(ω, Γ, Δh, Δg)

# Obtain data for ideal system
# Find the file "x.mat", "u.mat" and "xt.mat" from your corresponding file directory before executing following lines
vars=matread("x.mat")
vars1=matread("u.mat")
vars2=matread("xt.mat")
x=vars["x"]
x1_ideal=x[1,:];
x2_ideal=x[2,:];
u_ideal=vars1["u"];
xt=vars2["xt"]
x1t=xt[1,:];
x2t=xt[2,:];

l1_sol_x1=Vector{Float64}();
l1_sol_x2=Vector{Float64}();
append!(l1_sol_x1,-10);
append!(l1_sol_x2,10);
# L1 Reference sysem with perturbations AKA (non-implementable) ideal
# performance of the L1 system with full knowledge of the uncertainty
# ref_sys = reference_system(sys_p, ccm_p, ω, h, sigma, Sigma)
# ref_sol = solve(ref_sys, xs0, tspan, LambaEM(), progress = true, progress_steps = 1, saveat = 0.1, maxiters=5e7)

for i=1:20
    xs=[x1_ideal[i],x2_ideal[i]];
    us=u_ideal[i];
    ccm_p = ccm_params(xs, us, λ, W)
    # Sim time
    tspan=(0.01*i-0.01,0.01*i)
    # L1 + CCM system to handle disturbances
    l1_sys = l1_system(sys_p, ccm_p, l1_p, h, sigma, Sigma)
    l1_sol = solve(l1_sys, x0, tspan, LambaEM(), progress = true, progress_steps = 1, saveat = 1, maxiters=5e8)
    x0=[l1_sol[2][1],l1_sol[2][2]];
    append!(l1_sol_x1,l1_sol[2][1])
    append!(l1_sol_x2,l1_sol[2][2])

end
l1_sol_x1[21]=x1_ideal[21];

x0=[l1_sol_x1[21],l1_sol_x2[21]];
for j=21:180
    xs=[x1_ideal[j],x2_ideal[j]];
    us=u_ideal[j];
    ccm_p = ccm_params(xs, us, λ, W)
    # Sim time
    tspan=(0.01*j-0.01,0.01*j)
    # L1 + CCM system to handle disturbances
    l1_sys = l1_system(sys_p, ccm_p, l1_p, h, sigma, Sigma)
    l1_sol = solve(l1_sys, x0, tspan, LambaEM(), progress = true, progress_steps = 1, saveat = 1, maxiters=5e8)
    x0=[l1_sol[2][1],l1_sol[2][2]];
    append!(l1_sol_x1,l1_sol[2][1])
    append!(l1_sol_x2,l1_sol[2][2])

end

# MSE between actual system and reference system
mseactualandref=(l1_sol_x1-ref_sol_x1).^2+(l1_sol_x2-ref_sol_x2).^2;

# MSE between reference system and ideal system
mserefandideal=(ref_sol_x1).^2+(ref_sol_x2).^2;

# MSE between actual system and ideal system
mseactualandideal=(l1_sol_x1-x1_ideal).^2+(l1_sol_x2-x2_ideal).^2;

# Maximum 2-Wasserstein Distance between actual and ideal system
Wass2actualandideal=mseactualandideal.^0.5;
sort!(Wass2actualandideal,rev=true);
Wass2bound=Wass2actualandideal[1]

# Sample Ideal Trajectories
N=180;
xs1_sample=zeros(N+1,1000);
xs2_sample=zeros(N+1,1000);
xs1_sample[1,:].=-10;
xs2_sample[1,:].=10;
dt=0.01;
# 0.1 is the sqrt(dt)
for i= 1:N
    xs1_sample[i+1,:]=rand(Normal(x1_ideal[i+1],0.5*0.1*0.01*(x2_ideal[i]^2)),1000);
    xs2_sample[i+1,:]=rand(Normal(x2_ideal[i+1],2.0*0.1*0.01*(x2_ideal[i]^2)),1000);
end

# Sample True Trajectories
xt1_sample=zeros(N+1,1000);
xt2_sample=zeros(N+1,1000);
xt1_sample[1,:].=-10;
xt2_sample[1,:].=10;
for i= 1:N
    xt1_sample[i+1,:]=rand(Normal(l1_sol_x1[i+1],0.5*0.1*0.01*(l1_sol_x1[i]^2+l1_sol_x2[i]^2)),1000);
    xt2_sample[i+1,:]=rand(Normal(l1_sol_x2[i+1],2.0*0.1*0.01*(l1_sol_x1[i]^2+l1_sol_x2[i]^2)),1000);
end

#Generate Plot
p1=plot(l1_sol_x1,l1_sol_x2,label="CCM+L1+noise: Sample 1", title="Trajectory Plot: ω=50000, Γ=4e9",xlabel = L"x_1",ylabel = L"x_2")
plot!(x1_ideal,x2_ideal,label="Ideal+noise: Sample 1")
# define a function that returns a Plots.Shape
rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
plot!(rectangle(5,5.8,-8,0), opacity=.25,label="Obstacle Region")

# Plot 2000 Trajectories from ideal and actual system
p2=plot(xt1_sample[:,1],xt2_sample[:,1],label=false,linealpha=0.2,linewidth=0.1,title="2000 Total Sample True & Ideal Trajectories",xlabel = L"x_1",ylabel = L"x_2")
for k=2:1000
    plot!(xt1_sample[:,k],xt2_sample[:,k],label=false,linealpha=0.2,linewidth=0.1)
end
for l=1:1000
    plot!(xs1_sample[:,l],xs2_sample[:,l],label=false,linealpha=0.1,linewidth=0.1)
end
rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
plot!(rectangle(5,5.8,-8,0), opacity=.25,label="Obstacle Region")

# Plot the true sample trajectory using the input from ideal system
p3=plot(x1t,x2t,label="True sample trajectory using ideal input",xlabel = L"x_1",ylabel = L"x_2")
plot!(x1_ideal,x2_ideal,label="Ideal sample trajectory")
rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
plot!(rectangle(5,5.8,-8,0), opacity=.25,label="Obstacle Region")

# Save the graphs
savefig(p1,"Sample Trajectory Plot with Avoidance (omega = 50000, Gamma = 4e9).pdf")
savefig(p2,"2000 Total Sample ideal and actual trajectories (omega = 50000, Gamma = 4e9) .pdf")
savefig(p3,"Sample Trajectory Plot without l1-control (omega = 50000, Gamma = 4e9) .pdf")
