module TestParticleLimit

export timesample, ionize_stepped
export N0_kinetx, dσ_kinetx, z_R1, z_R2, v_kinetx, vx_kinetx
export N⃗_excitation, k_excitation, k_equilibrium
export Nn_equilibrium, dNn_equilibrium, Nn_excitation, dNn_excitation
export Ni_equilibrium, dNi_equilibrium, Ni_excitation, dNi_excitation
export gaussian_sample, moveion, crres_sample
export propagate_ions, final_ions
export N_madeup, k_madeup, dNi_madeup

using Unitful
using StaticArrays
using LinearAlgebra

include("TruncatedNormals.jl")

using .TruncatedNormals
export TruncatedNormal

"""
    timesample(τ)
    timesample(τ1, τs...)

Produce a random sample of photoionization time. More precicely, sample from a
[hypoexponential distribution](https://en.wikipedia.org/wiki/Hypoexponential_distribution)
with rate parameters λ_i = 1/τ[i].

In a photoequlibrium model the time until a particular neutral is ionized is memoryless in
the sense that the expected time until photoionization does not depend on how long the
neutral has already been exposed to the light. Therefore, photoionization time has the
exponetial distribution. (Proofs are readily available online.)

We also allow ionization to occur in a sequence of memoryless stages. In particular the
model we're using for gradual distriubuted ionization is that X converts to neutral barium
with a time constant of about 5 seconds and neutral barium then photoionizes with a time
constant of about 30 seconds. So the total time to ionization is distributed as a
hypoexponential and we can sample from that distribution with `timesample(5u"s", 30u"s")`.

It's possible this (or a similar phase-type distribution/markov process) could be used to
model the process of electron energization through several metastable states until
photoionization occurs.
"""
function timesample(τ1, τs...)
    sum(τ * randexp() for τ in Iterators.flatten((τ1, τs)))
end

const k1i = 0.0014u"s^-1"
const k2i = 0.033u"s^-1"
const kti = 0.071u"s^-1"
const k12 = 0.58u"s^-1"
const k21 = 4.05u"s^-1"
const kt1 = 0.4u"s^-1"
const k1t = 0.55u"s^-1"
const k2t = 12.19u"s^-1"
const kt2 = 2.48u"s^-1"
const Q = [ # Transition Matrix (transposed)
    0u"s^-1"       k1i            k2i            kti;
    0u"s^-1" -(k1i+k12+k1t)       k21            kt1;
    0u"s^-1"       k12      -(k2i+k21+k2t)       kt2;
    0u"s^-1"       k1t            k2t      -(kti+kt1+kt2)
]

const τ_equlibrium = 28u"s"
const N0_kinetx = 2.631e24
const FWHMrate_kinetx = 3u"km/s"
const dσ_kinetx = FWHMrate_kinetx/2√(2log(2))
const v_kinetx = 2.5u"km/s"
const vx_kinetx = 2.153u"km/s"
const v_crres = 9.5u"km/s"
const v0_crres = 1.33u"km/s"
const vth_crres = 0.296u"km/s"
const dist_crres = TruncatedNormal(v0_crres, vth_crres*√2, zero(v0_crres))
const z_R1 = 4.13u"km"
const z_R2 = 2.81u"km"

# Made up model
# the parameter α adjusts the excitaiton rate. The long time ionization rate
# is also affected, but only slightly. I'm not sure how to correct for that.
function Q_madeup(α)
    _k1i = k1i
    _k2i = k2i
    _kti = kti
    _k12 = k12*α
    _k21 = k21*α
    _kt1 = kt1*α
    _k1t = k1t*α
    _k2t = k2t*α
    _kt2 = kt2*α
    [ # Transition Matrix (transposed)
        0u"s^-1"       _k1i            _k2i            _kti;
        0u"s^-1" -(_k1i+_k12+_k1t)       _k21            _kt1;
        0u"s^-1"       _k12      -(_k2i+_k21+_k2t)       _kt2;
        0u"s^-1"       _k1t            _k2t      -(_kti+_kt1+_kt2)
    ]
end
N⃗_madeup(t, α::Number; N0=N0_kinetx) = N⃗_madeup(t, Q_madeup(α); N0)
function N⃗_madeup(t, Q::AbstractMatrix; N0=N0_kinetx)
    exp(ustrip.(Unitful.NoUnits, t*Q)) * [0,N0,0,0]
end

function k_madeup(t, α)
    N⃗ = N⃗_madeup(t, α)
    Nn = @views sum(N⃗[2:end])
    α = N⃗[2]/Nn
    β = N⃗[3]/Nn
    γ = N⃗[4]/Nn
    k1i*α + k2i*β + kti*γ
end
dNi_madeup(t, α::Number; N0=N0_kinetx) = dNi_madeup(t, Q_madeup(α); N0)
function dNi_madeup(t, Q::AbstractMatrix; N0=N0_kinetx)
    N0 * (Q * exp(t*Q))[1,2]
end

Nn_equilibrium(t; N0=N0_kinetx) = N0 * exp(-t/τ_equlibrium)
dNn_equilibrium(t; N0=N0_kinetx) = -N0/τ_equlibrium * exp(-t/τ_equlibrium)
k_equilibrium(t) = 1/τ_equlibrium

N⃗_excitation(t; N⃗0=SA[0,N0_kinetx,0,0]) = exp(ustrip.(Unitful.NoUnits, t*Q)) * N⃗0
function k_excitation(t)
    N⃗ = N⃗_excitation(t)
    Nn = @views sum(N⃗[2:end])
    α = N⃗[2]/Nn
    β = N⃗[3]/Nn
    γ = N⃗[4]/Nn
    k1i*α + k2i*β + kti*γ
end
Nn_excitation(t; N0=N0_kinetx) = @views N0 * sum(exp(t*Q)[2:end,2])
dNn_excitation(t; N0=N0_kinetx) = @views N0 * sum((Q * exp(t*Q))[2:end,2])

Ni_equilibrium(t; N0=N0_kinetx) = N0 * (1 - exp(-t/τ_equlibrium))
dNi_equilibrium(t; N0=N0_kinetx) = N0/τ_equlibrium * exp(-t/τ_equlibrium)

Ni_excitation(t; N0=N0_kinetx) = N0 * exp(t*Q)[1,2]
dNi_excitation(t; N0=N0_kinetx) = N0 * (Q * exp(t*Q))[1,2]

function gaussian_sample(t, n; dσ=dσ_kinetx, v=vx_kinetx)
    σ = t*dσ
    x0 = v*t
    T = SVector{3, promote_type(typeof(σ), typeof(x0))}
    result = Vector{T}(undef, n)
    for i in eachindex(result)
        result[i] = T(σ*randn() + x0, σ*randn(), σ*randn())
    end
    result
end

function sphere_sample()
    normalize(SA[randn(), randn(), randn()])
end

function crres_sample(t, n; v=v_crres, v0=v0_crres, vth=vth_crres)
    dist = TruncatedNormal(v0, vth*√2, zero(v0))
    x0 = v*t
    T = SVector{3, promote_type(eltype(dist), typeof(x0))}
    result = Vector{T}(undef, n)
    for i in eachindex(result)
        d = sphere_sample()
        vv = rand(dist)
        result[i] = T(t*vv*d[1]+x0, t*vv*d[2], t*vv*d[3])
    end
    result
end

stepcenters(stepedges::AbstractRange) = range(step(stepedges)/2, last(stepedges)+step(stepedges)/2, length(stepedges))

"""
    ionize_stepped(times, dNi, α, sample_ions)

At each time in times, produce a list of macroparticle ions born in each timestep.

α is the number of macroparicles for N0 microparticles where N0 is the initial number of neutrals.
I.e. β = α/N0. I.e. α is roughly the total number of particles (i.e. `sum(length, ionize_stepped(...))` )
in the limit of long duration and short timestep.
"""
function ionize_stepped(times::AbstractRange, dNi, α, sample_ions, N0=N0_kinetx)
    dt = step(times)
    ion_times = stepcenters(times) # dt == step(times) ≈ step(centers) the discretized time of ionization
    num_new_ions = (round(Int, α/N0*dt*dNi(t; N0)) for t in ion_times)
    [sample_ions(t, n) for (t,n) in zip(ion_times, num_new_ions)]
end

function moveion(ion, t0, tf)
    if iszero(t0)
        error("We don't deal with ions \"born\" at t=0")
    end

    _moveion(ion, t0, tf)
end

function _moveion(ion::StaticVector{3}, t0, tf)
    T = similar_type(ion)
    T(
        ion[1],
        ion[2],
        ion[3] * tf/t0
    )
end
function _moveion(ion::AbstractVector, t0, tf)
    result = similar(ion)
    result[1] = ion[1]
    result[2] = ion[2]
    result[3] = ion[3] * tf/t0
    result
end

"""
    propagate_ions(tfinal, times, ions)

Calculates the positions of the ions at the time `tfinal` in the small gyroradius and test particle limits.
The `ions` argument should be the ions born at the times in `times`. I.e. ions = ionize_stepped(times, ...).
See final_ions() for a typical use case.
"""
function propagate_ions(tfinal, times::AbstractRange, ions)
    ion_times = stepcenters(times)
    first(ion_times)-step(ion_times)/2 <= tfinal <= last(ion_times)+step(ion_times)/2 || error("tfinal must be within the range of times")

    tfinal_index = searchsortednearest(ion_times, tfinal)
    result_length = @views sum(length, ions[1:tfinal_index])
    result = Vector{eltype(eltype(ions))}(undef, result_length)

    i = 1
    for (t0, g) in @views zip(ion_times[1:tfinal_index], ions[1:tfinal_index])
        for ion in g
            result[i] = moveion(ion, t0, tfinal)
            i += 1
        end
    end
    if i != length(result) + 1
        error("Unexpected Error: result hasn't been filled correctly.")
    end
    result
end

function final_ions(tfinal, nsteps, ionizationrate, cloudsampler, macroparticles)
    times = range(zero(tfinal), tfinal, length=nsteps)
    ions = ionize_stepped(times, ionizationrate, macroparticles, cloudsampler)
    propagate_ions(tfinal, times, ions)
end

function searchsortednearest(a, x)
    idx = searchsortedfirst(a, x)
    if idx==1
        return idx
    end

    if idx > length(a)
        return length(a)
    end

    if a[idx] == x
        return idx
    end

    if abs(a[idx]-x) < abs(a[idx-1]-x)
        return idx
    else
        return idx-1
    end
end

include("analytic.jl")

end # module TestParticleLimit
