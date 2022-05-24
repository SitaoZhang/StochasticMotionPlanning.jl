using DifferentialEquations
struct SysParams{fF,BF}
    f::fF
    B::BF
end
sys_params(f::Function, B::Function) = SysParams(f, B)
function sys_params(f::Function, B::Union{Matrix,Vector,UniformScaling})
    Bfcn(x) = B
    sys_params(f, Bfcn)
end
function sys_params(f::Vector, B)
    Ffcn(x) = f
    sys_params(Ffcn, B)
end

struct CCMParams{T1,xsF,usF,WF}
    xs::xsF
    us::usF
    λ::T1
    W::WF
    A::Matrix{T1}
    T::Matrix{T1}
    Ts::Matrix{T1}
    n::Vector{T1}
    w::Vector{T1}
    γ0::Matrix{T1}
end
function ccm_params(xs::Function, us::Function, λ, W::Function; D=5, N=7)
    nodes, weights = _chebpts(N)
    T = _chebpoly(nodes, D)
    Ts = _chebpolyder(T, nodes, D)
    Ti = T[:,2:N-1]; Te = T[:,[1,N]]
    A = inv([2*Ti*Ti' Te;Te' zeros(2,2)])
    γ0 = zeros(length(xs(0)), N-2)
    λ = convert(eltype(T),λ)
    CCMParams(xs, us, λ, W, A, T, Ts, nodes, weights, γ0)
end
function ccm_params(xs::Vector, us, λ, W; kwargs...)
    xf(t) = xs
    ccm_params(xf, us, λ, W; kwargs...)
end
function ccm_params(xs::Function, us::Union{Vector,Real}, λ, W; kwargs...)
    uf(t) = us
    ccm_params(xs, uf, λ, W; kwargs...)
end


struct FlatCCMParams{T,xsF,usF,WT}
    xs::xsF
    us::usF
    λ::T
    M::WT
end
function ccm_params(xs, us, λ, W)
    FlatCCMParams(xs, us, λ, inv(W))
end


struct L1Params{T,PT,AmT}
    ω::T
    Γ::T
    Δh::T
    Δg::T
    P::PT
    Am::AmT
end
function l1_params(ω, Γ, Δh,Δg, P = I, Am = -I)
    ω, Γ, Δh, Δg = promote(ω, Γ, Δh, Δg)
    L1Params(ω, Γ, Δh,Δg, P, Am)
end

W1=WienerProcess(0.0,2.0)
W2=WienerProcess(0.0,0.5)

struct NominalProblem{A}
    prob::SDEProblem{A}
end
function nominal_system(sys, ccm, Sigma)
    NominalProblem(SDEProblem(_nominal_system!, _nominal_noise!,nothing, nothing, (sys, ccm, Sigma);noise=W1))
end

function nominal_system(sys, ccm, h, Sigma)
    NominalProblem(SDEProblem(_perturbed_system!, _nominal_noise!, nothing, nothing, (sys, ccm, h, Sigma);noise=W1))
end


struct ReferenceProblem{A}
    prob::SDEProblem{A}
end
function reference_system(sys, ccm, ω, h, Sigma, sigma)
    ReferenceProblem(SDEProblem(_reference_system!, _reference_noise!, nothing, nothing, (sys, ccm, ω, h, Sigma, sigma);noise=W1))
end
reference_system(sys, ccm, l1::L1Params, h, Sigma, sigma) = reference_system(sys, ccm, l1.ω, h, Sigma, sigma)


struct L1Problem{A}
    prob::SDEProblem{A}
end
function l1_system(sys, ccm, l1, h, Sigma, sigma)
    L1Problem(SDEProblem(_l1_system!,_l1_noise!, nothing, nothing, (sys, ccm, l1, h, Sigma, sigma);noise=W2))
end
