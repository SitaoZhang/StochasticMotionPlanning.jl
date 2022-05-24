function _nominal_system!(D, z, (sys, ccm, Sigma), t)
    @unpack x = z
    @unpack f, B = sys

    uc = _u_ccm(x, sys, ccm, t)
    D.x = f(x) + B(x)*uc
end

function _nominal_noise!(D, z, (sys, ccm, Sigma), t)
    @unpack x = z
    @unpack f, B = sys

    D.x = B(x)*Sigma(t, x)
end

function _perturbed_system!(D, z, (sys, ccm, h), t)
    @unpack x = z
    @unpack f, B = sys

    uc = _u_ccm(x, sys, ccm, t)
    D.x = f(x) + B(x)*(uc + h(t, x))
end

function _reference_system!(D, z, (sys, ccm, ω, h, Sigma, sigma), t)
    @unpack x, η = z
    @unpack f, B = sys

    uc = _u_ccm(x, sys, ccm, t)
    u = uc - η

    D.η = -ω*η + ω*h(t, x)
    D.x = f(x) + B(x)*(u + h(t, x))
end

function _reference_noise!(D, z, (sys, ccm, ω, h, Sigma, sigma), t)
    @unpack x, η = z
    @unpack f, B = sys

    Sigma_r = Sigma(t,x) - η
    D.η = -ω*η + ω*sigma(t, x)
    D.x = B(x)*(Sigma_r + sigma(t,x))
end

function _l1_system!(D, z, (sys, ccm, l1, h, Sigma, sigma), t)
    @unpack x, xhat, σhat, σhat_1, ηhat, ηhat_1 = z
    @unpack f, B = sys

    uc = _u_ccm(x, sys, ccm, t)
    ua = _u_l1!(D, z, sys, l1, uc)
    u = uc + ua

    D.x = f(x) + B(x)*(u + h(t, x))
end

function _l1_noise!(D, z, (sys, ccm, l1, h, Sigma, sigma), t)
    @unpack x, xhat, σhat, σhat_1, ηhat, ηhat_1 = z
    @unpack f, B =sys

    D.x = B(x)*(Sigma(t,x)+sigma(t, x))
end

function _initialize_geodesic!(ccm::CCMParams, x0)
	γ0(s) = ccm.xs(0)*(1-s) + x0*s
	ccm.γ0 .= hcat(γ0.(ccm.n[2:end-1])...)
end

_initialize_geodesic!(ccm::FlatCCMParams, x0) = nothing

function solve(sys::NominalProblem, x0, tspan, varargs...; kwargs...)
    _initialize_geodesic!(sys.prob.p[2], x0)
    z0 = ComponentArray(x=x0)
    solve(sys.prob, varargs...; u0 = z0, tspan = tspan, kwargs...)
end

function solve(sys::ReferenceProblem, x0, tspan, varargs...; kwargs...)
    _initialize_geodesic!(sys.prob.p[2], x0)
    η0 = sys.prob.p[4](0, x0)
    z0 = ComponentArray(x=x0, η=η0)
    solve(sys.prob, varargs...; u0 = z0, tspan = tspan, kwargs...)
end

function solve(sys::L1Problem, x0, tspan, varargs...; kwargs...)
    _initialize_geodesic!(sys.prob.p[2], x0)
    σ0 = sys.prob.p[4](0, zero(x0))
    z0 = ComponentArray(x=x0, xhat=x0, σhat=σ0, σhat_1=σ0, ηhat=σ0, ηhat_1=σ0)
    solve(sys.prob, varargs...; u0 = z0, tspan = tspan, kwargs...)
end
