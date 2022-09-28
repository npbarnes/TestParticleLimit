module TruncatedNormals
using SpecialFunctions
using Random

export TruncatedNormal

Φ(x) = 1/2 * (1 + erf(x/√2))
Φinv(p) = √2 * erfinv(2p-1)

struct TruncatedNormal{T}
    μ::T
    σ::T
    a::T
end
TruncatedNormal(a,b,c) = TruncatedNormal(promote(a,b,c)...)

struct Regime1 end
struct Regime2 end
struct Regime3 end

Base.eltype(::Type{TruncatedNormal{T}}) where T = T

function Random.Sampler(::Type{<:AbstractRNG}, d::TruncatedNormal, ::Random.Repetition)
    ā = (d.a - d.μ)/d.σ
    if ā < -4
        return Random.SamplerSimple(d, (Regime1(), d.μ, d.σ, ā))
    elseif ā > 4
        return Random.SamplerSimple(d, (Regime3(), d.μ, d.σ, ā, ā^2))
    else
        t = Φ(ā)
        tc = 1-t
        return Random.SamplerSimple(d, (Regime2(), d.μ, d.σ, t, tc))
    end
end

Random.rand(rng::AbstractRNG, sp::Random.SamplerSimple{<:TruncatedNormal}) = sample_truncnorm(rng, sp.data...)

function sample_truncnorm(rng, ::Regime1, μ, σ, ā)
    x̄ = randn(rng)
    while x̄ < ā
        x̄ = randn(rng)
    end
    σ*x̄ + μ
end

function sample_truncnorm(rng, ::Regime2, μ, σ, t, tc)
    u = tc*rand(rng) + t
    x̄ = Φinv(u)
    σ*x̄ + μ
end

function sample_truncnorm(rng, ::Regime3, μ, σ, ā, ā2)
    u = rand(rng)
    x̄ = √(ā2 - 2log(1-u))
    v = rand(rng)
    while v <= x̄/ā
        u = rand(rng)
        x̄ = √(ā2 - 2log(1-u))
        v = rand(rng)
    end
    σ*x̄ + μ
end

end # module
