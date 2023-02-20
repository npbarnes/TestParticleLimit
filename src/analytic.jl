
using QuadGK
using UnPack
using PhysicalConstants.CODATA2018: m_p


export m_Ba
export P_thermal, nn_payload, 𝓅_payload, ni_observed

const m_Ba = 137m_p

P_thermal(v; dσ=dσ_kinetx)  = 4π*(1/(2π*dσ^2))^(3/2) * v^2 * exp(-(1/2)*(v/dσ)^2)

function nn_payload(t,r; parameters)
    if r<0u"km"
        error("Unexpected negative radius")
    end
    if t ≈ zero(t) || t < zero(t)
        return 0.0u"km^-3"
    end
    @unpack Nn, N0, P, dσ = parameters
    Nn(t; N0)*P(r/t; dσ)/(4π*r^2*t)
end

"""
    𝓅(t,r; parameters)

ion production rate density. I.e. 𝓅 = ∂nᵢ/∂t + ∇⋅(nᵢuᵢ) = k(t) * nₙ
"""
function 𝓅_payload(t,r; parameters)
    @unpack k = parameters
    k(t) * nn_payload(t,r; parameters)
end
𝓅_payload(t,x,z; parameters) = 𝓅_payload(t, √(x^2 + z^2); parameters)

function ni_observed(t; parameters)
    @unpack vD, zR = parameters
    xi(s) = vD*(t-s)
    zi(s) = zR*s/t
    integrand(s) = 𝓅_payload(s, xi(s), zi(s); parameters)
    result, err = quadgk(integrand, zero(t), t)
    result
end

function f_1d_drift_observed(t, vp′; parameters)
    @unpack vD, zR = parameters
    if t ≈ zero(t) || t < zero(t) || vp′ < vD
        return 0.0u"s * km^-4"
    end
    t*vD/vp′^2 * 𝓅_payload(t*vD/vp′, t*vD*(1 - vD/vp′), zR*vD/vp′; parameters)
end

function f_2d_drift_observed(t, v⃗p′; parameters)
    vx′, vy′ = v⃗p′
    vp′ = √(vx′^2 + vy′^2)
    f_1d_drift_observed(t, vp′)/(2π*vp′)

end

export f_2d
"""
Phase space density integrated over all z. The z dependence is a delta function
"""
function f_2d(t, v⃗perp; parameters)
    @unpack vD, zR = parameters
    vx, vy = v⃗perp
    vp = √((vx + vD)^2 + vy^2)
    result = t*vD/(2π*vp^3) * 𝓅_payload(t*vD/vp, t*vD*(1-vD/vp), zR*vD/vp; parameters)
    vp < vD ? zero(result) : result
end

export f_E
function f_E(t, E, θ; parameters)
    v = √(2E/m_Ba) |> u"km/s"
    f_2d(t, (-v*cos(θ), v*sin(θ)); parameters)
end

function differential_intensity(t, v⃗; parameters)
    @unpack zR = parameters
    vx, vy = v⃗
    vp = √(vx^2 + vy^2)
    vz = zR/t
    v = √(vx^2 + vy^2 + vz^2)

    f = f_2d(t, v⃗; parameters)
    I = f*v^2/m_Ba

end

function _F(t; parameters)
    t_I(t,x; vD) = t - x/vD
    z_I(t,x; zR, vD) = zR - zR*x/(t*vD)
    r_I(t,x; zR, vD) = √(x^2 + z_I(t,x; zR, vD)^2)
    @unpack P, Nn, vD, zR = parameters
    # Parameterize the ionization manifold as γ(s) = [s, 0, zR - (zR*s)/(t*vD)].
    # Then J = ||γ'(s)||. F(t) = ∫dni(t_I(s),γ(s)) ||γ'(s)|| ds
    J = √(1 + (zR/(t*vD))^2)
    integrand(s) = J * 𝓅_payload(t_I(t,s; vD), r_I(t,s; zR, vD); parameters)
    result, err = quadgk(integrand, 0u"km", t*vD)
    result
end

export _ni
function _ni(t; parameters)
    @unpack zR, vD = parameters
    _F(t; parameters)/√((zR/t)^2 + vD^2)
end

export production_vs_vp
function production_vs_vp(t, vp; parameters)
    @unpack zR, vD = parameters
    xx = vD*t*(vp - vD)/vp
    if xx < 0u"km" || xx > vD*t
        return 0u"km^-3 * s^-1"
    end
    𝓅_payload(t_I(t,xx; vD), r_I(t,xx; zR, vD); parameters)
end

export direct_ion_density
function direct_ion_density(t; parameters)
    @unpack vD, zR = parameters
    # parameterize ionization location in terms of ionization time
    xi(s) = vD*(t-s)
    zi(s) = zR*s/t
    ri(s) = √(xi(s)^2 + zi(s)^2)
    integrand(s) = 𝓅_payload(s, ri(s); parameters)
    result, err = quadgk(integrand, zero(t), t)
    result
end

parameters = (
    # Base parameters (values)
    N0=N0_kinetx,
    dσ=dσ_kinetx,
    zR=z_R2,
    vD=v_kinetx,

    # Function parameters
    P=P_thermal,
    Nn=Nn_excitation,
    k=k_excitation,
)
