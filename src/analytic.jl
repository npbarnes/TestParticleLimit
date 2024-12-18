
using Cubature
using QuadGK
using UnPack
using PhysicalConstants.CODATA2018: m_p

export m_Ba
export P3_thermal, P_thermal, nn_payload, 𝓅_payload, ni_observed

const m_Ba = 137m_p

struct _Zero end
_unitless_abstol(::_Zero, u) = 0
_unitless_abstol(a, u) = ustrip(u, a)

"""
The Cubatrue package does not support units, so this function is a workaround that strips
units first, then does the dimensionless integration, and reattaches appropriate units
before returning the result. This implementation has some limitations:
1) The integrand f is called once on an arbitrary point inside the integration domain
in order to determine the units. This may have implications for functions with side effects,
or computationally expensive functions, or functions that happen to have a sigularity at the
representative point. Those types of functions most likely would be problematic for cubature
anyway.
2) If the unit of the value returned by f(x) depends on the value of x (not just the units of x),
then ucubature will produce the correct result with units, but it might not be the units you expect.
The result of f(x) will be converted to unit(f(x_rep)) before being stripped.
"""
uhcubature(args...; kwargs...) = _ucubature(hcubature, args...; kwargs...)
upcubature(args...; kwargs...) = _ucubature(pcubature, args...; kwargs...)
function _ucubature(cubature::Function, f::Function, xmin, xmax; reltol=1e-8, abstol=_Zero(), maxevals=0)
    x = (xmin .+ xmax) ./ 2
    # We implicitly convert to preferred units when we compute the representative point, x.
    xunits = unit.(x)
    funit = unit(f(x))
    Iunit = funit * prod(xunits)

    clean_abstol = _unitless_abstol(abstol, Iunit)
    clean_f(args) = ustrip(funit, f(args .* xunits))
    clean_xmin = ustrip.(xmin)
    clean_xmax = ustrip.(xmax)
    clean_int, clean_err = cubature(clean_f, clean_xmin, clean_xmax; reltol, abstol=clean_abstol, maxevals)
    (clean_int*Iunit, clean_err*Iunit)
end

P3_thermal(v::AbstractArray; dσ=dσ_kinetx) = P3_thermal(norm(v); dσ)
function P3_thermal(v::Number; dσ=dσ_kinetx)
    _dσ = uconvert(unit(v), dσ)
    (1/2π)^(3/2) * _dσ^-3 * exp(-(1/2)*(v/dσ)^2)
end

function P_thermal(v; dσ=dσ_kinetx)
    result = 4π * v^2 * P3_thermal(v; dσ)

    # In the limit as v -> inf, P -> 0
    isinf(v) ? zero(result) : result
end

export nn_payload

"""
    nn_payload(t,r; parameters)

Neutral density as a function of time and radius (in the payload frame). Gives the
correct values when either  t is zero, or r in zero, or both are zero. Undefined when
either t or r is negative.
"""
function nn_payload(t,r; parameters)
    if r == zero(r) && t == zero(t)
        return Inf*unit(r)^-3
    elseif t == zero(t)
        return zero(r^-3)
    end

    @unpack Nn, N0, P3, dσ = parameters
    Nn(t; N0) * t^-3 * P3(r/t; dσ)
end

"""
    𝓅(t,r; parameters)

ion production rate density. I.e. 𝓅 = ∂nᵢ/∂t + ∇⋅(nᵢuᵢ) = k(t) * nₙ
"""
function 𝓅_payload(t,r; parameters)
    @unpack k = parameters
    k(t) * nn_payload(t,r; parameters)
end
𝓅_payload(t,x ,y, z; parameters) = 𝓅_payload(t, √(x^2 + y^2 + z^2); parameters)

export 𝓅_ambient
function 𝓅_ambient(t, x′, y′, z′; parameters)
    @unpack vD = parameters
    𝓅_payload(t, x′ - t*vD, y′, z′; parameters)
end

export density1_ambient
function density1_ambient(t, x′, y′, z′; parameters)
    integrand(tᵢ) = tᵢ/t * 𝓅_ambient(tᵢ, x′, y′, z′*tᵢ/t; parameters)
    quadgk(integrand, zero(t), t)
end

export density1_payload
function density1_payload(t, x, y, z; parameters)
    @unpack vD = parameters
    density1_ambient(t, x+t*vD, y, z; parameters)
end

export f_2d_ambient
function f_2d_ambient(t, x′, y′, z′, vx′, vy′; parameters)
    ρ = √(x′^2 + y′^2)
    vp = √(vx′^2 + vy′^2)

    # Handle edge cases:
    if t<zero(t) || isapprox(t, zero(t); atol=eps(1.0u"s"))
        return 0.0u"km^-5*s^2"
    end
    if vp <= ρ/t
        return 0.0u"km^-5*s^2"
    end

    ρ / (2π*vp^3) * 𝓅_ambient(ρ/vp, x′, y′, z′*ρ/(vp*t); parameters)
end

export f_2d_payload
function f_2d_payload(t, x, y, z, vx, vy; parameters)
    @unpack vD = parameters
    f_2d_ambient(t, x + t*vD, y, z, vx + vD, vy; parameters)
end

export f_2d_payload_cylindrical
function f_2d_payload_cylindrical(t, x, y, z, vxy, vθ; parameters)
    #vz = z/t
    #v = √(vxy^2 + vz^2)
    vx = vxy*cos(vθ)
    vy = vxy*sin(vθ)
    f_2d_payload(t, x, y, z, vx, vy; parameters)
end

export density2_payload
function density2_payload(t, x, y, z; parameters)
    function integrand(args)
        vx = args[1]
        vy = args[2]

        #vz = z/t
        #v = √(vxy^2 + vz^2)

        f_2d_payload(t,x,y,z, vx, vy; parameters)
    end
    uhcubature(integrand, [-10, -10]u"km/s", [10,10]u"km/s"; reltol=1e-3)
end

export density3_payload
function density3_payload(t, x, y, z; parameters)
    function integrand(args)
        vxy = args[1]
        vθ = args[2]
        vxy * f_2d_payload_cylindrical(t, x, y, z, vxy, vθ; parameters)
    end
    uhcubature(integrand, [0u"km/s", deg2rad(-180)], [10u"km/s", deg2rad(180)]; reltol=1e-3)
end

#= using PyPlot: plt
export plot_dist_xy
function plot_dist_xy(parameters)
    fig, ax = plt.subplots()
    vxs = range(-10u"km/s", 10u"km/s", length=99)
    vys = range(-10u"km/s", 10u"km/s", length=99)

    @unpack vD = parameters
    t = 1u"s"
    Δ = vD*t
    fs = [f_2d_ambient(1u"s", -Δ, 0u"km", 3u"km", vx, vy; parameters) for vx in vxs, vy in vys]

    ax.pcolormesh(ustrip.(u"km/s", vxs), ustrip.(u"km/s", vys), transpose(ustrip.(u"s^2 * km^-5", fs)))
    ax.set_aspect("equal")
    fig, ax
end

export plot_dist_cylindrical
function plot_dist_cylindrical(parameters)
    fig, ax = plt.subplots(subplot_kw=Dict("projection"=>"polar"))
    vs = range(0u"km/s", 10u"km/s", length=100)
    vts = range(-π, π, length=99)

    fs = [f_2d_payload_cylindrical(1u"s", 0u"km", 2u"km", 3u"km", v, vt; parameters) for vt in vts, v in vs]

    ax.pcolormesh(vts, ustrip.(u"km/s", vs), transpose(ustrip.(u"s^2 * km^-5", fs)))

    fig, ax
end =#

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
