using UnPack
using Unitful

const Ω = 1.0
const N0 = 1.0
const σx = 1.0
const σv = 1.0

struct R
    x0::Float64
    y0::Float64
    z0::Float64
    vx0::Float64
    vy0::Float64
    vz0::Float64
end
R(x,v) = R(x..., v...)

struct Ξ
    x::Float64
    y::Float64
    z::Float64
    vx::Float64
    vy::Float64
    vz::Float64
end
Ξ(x,v) = Ξ(x..., v...)

position(ξ::Ξ) = SA[ξ.x, ξ.y, ξ.z]
position(r::R) = SA[r.x0, r.y0, r.z0]

velocity(ξ::Ξ) = SA[ξ.vx, ξ.vy, ξ.vz]
velocity(r::R) = SA[r.vx0, r.vy0, r.vz0]

function g3(ζ, σ)
    (2π*σ)^(-3/2) * exp(-1/2 * (ζ⋅ζ/σ^2))
end

function x(r,s)
    @unpack vy0, vx0, x0 = r
    -vy0/Ω * cos(Ω*s) + vx0/Ω * sin(Ω*s) + vy0/Ω + x0
end

function y(r,s)
    @unpack vx0, vy0, vx0, y0 = r
    vx0/Ω * cos(Ω*s) + vy0/Ω * sin(Ω*s) - vx0/Ω + y0
end

function z(r,s)
    @unpack z0,vz0 = r
    z0 + vz0*s
end

function vx(r,s)
    @unpack vy0,vx0 = r
    vy0*sin(Ω*s) + vx0*cos(Ω*s)
end

function vy(r,s)
    @unpack vx0, vy0 = r
    -vx0*sin(Ω*s) + vy0*cos(Ω*s)
end

function vz(r,s)
    r.vz0
end

function ξ(r,s)
    Ξ(
        x(r,s),
        y(r,s),
        z(r,s),
        vx(r,s),
        vy(r,s),
        vz(r,s)
    )
end

function x0(ξ, t)
    @unpack x,vx,vy = ξ
    x - vx/Ω * sin(Ω*t) + vy/Ω * (1-cos(Ω*t))
end

function y0(ξ, t)
    @unpack y, vx, vy = ξ
    y + vx/Ω * (cos(Ω*t)-1) - vy/Ω * sin(Ω*t)
end

function z0(ξ, t)
    @unpack z, vz = ξ
    z - vz*t
end

function vx0(ξ, t)
    @unpack vx,vy = ξ
    vx*cos(Ω*t) - vy*sin(Ω*t)
end

function vy0(ξ, t)
    @unpack vx,vy = ξ
    vx*sin(Ω*t) + vy*cos(Ω*t)
end

function vz0(ξ, t)
    ξ.vz
end

function r(ξ, t)
    R(
        x0(ξ,t),
        y0(ξ,t),
        z0(ξ,t),
        vx0(ξ,t),
        vy0(ξ,t),
        vz0(ξ,t)
    )
end

struct Gaussian
    μ::SVector{3,Float64}
    σ::Float64
end
function (g::Gaussian)(q)
    @unpack μ, σ = g
    (1/(2π*σ^2))^(3/2) * exp(-1/2 * (x-μ)⋅(x-μ)/σ^2)
end

P6(ξ::Ξ) = P6(position(ξ), velocity(ξ))
function P6(x,v; v_rocket)
    x_gauss = Gaussian(zeros(3), 0.1)
    v_gauss = Gaussian(SA[v_rocket, 0.0, 0.0], dσ_kinetx)

    x_gauss(x) * v_gauss(v)
end

function fn(t,ξ; Nn)
    ξ2 = Ξ(position(ξ) - t*velocity(ξ), velocity(ξ))
    Nn(t) * P6(ξ2)
end

function h(r,s)

end

fi(t,ξ) = h(r(ξ,t), t)
