function _proj(Δθ, θ, y; ϵ = 0.1)
    f = ((ϵ + 1)*θ'*θ - Δθ^2)/(ϵ*Δθ^2)
    ∇f = (2*(ϵ + 1)/(ϵ*Δθ^2))*θ

    if f ≥ 0 && ∇f'*y > 0
        return y - ∇f*(∇f'*y)*f/(∇f'*∇f)
    else
        return y
    end
end

function _u_l1!(D, z, sys, l1, uc)
    @unpack x, xhat, σhat, σhat_1, ηhat, ηhat_1 = z
    @unpack f, B = sys
    @unpack ω, Γ, Δh, Δg, P, Am = l1

    ua = -ηhat
    xtilde = xhat - x

    D.xhat = f(x) + B(x)*(uc + ua + σhat) + Am*xtilde
    D.σhat = Γ* _proj(Δh, σhat, -B(x)'*P*xtilde)
    D.σhat_1 = Γ* _proj(Δg, σhat_1, -B(x)'*P*xtilde)
    D.ηhat = -ω*ηhat + ω*σhat
    D.ηhat_1 = -ω*ηhat_1 + ω*σhat_1
    return ua
end
