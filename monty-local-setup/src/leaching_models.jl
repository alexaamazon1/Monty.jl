module LeachingModels

using Random: AbstractRNG, Xoshiro
using UnPack: @unpack
using Unitful: Quantity, ustrip, @u_str, unit

using DocStringExtensions

#------------------------------------------------------------------------------
export AbstractLeachingModel

abstract type AbstractLeachingModel{𝒯} end

function (L::AbstractLeachingModel{𝒯})(t::Quantity{𝒯}) where {𝒯}
    @assert unit(t) == u"yr"
    L(t |> ustrip)
end

#------------------------------------------------------------------------------
export NoLeaching
export ExponentialLeaching, MultiExponentialLeaching, SeasonalLeaching

#--------------------------------------

function leachnoise(rng::AbstractRNG, 𝓁::𝒯, σ::𝒯)::𝒯 where {𝒯}
    iszero(σ) && return 𝓁
    x = log(𝓁 / (1 - 𝓁))
    x += σ * randn(rng)
    1 / (1 + exp(-x))
end

#--------------------------------------

"""Leaching model with zero leaching at any time"""
struct NoLeaching{𝒯} <: AbstractLeachingModel{𝒯}
    function NoLeaching()
        new{Float64}()
    end
end

NoLeaching(T::Type) = NoLeaching{T}()

function (w::NoLeaching{𝒯})(t::𝒯)::𝒯 where {𝒯}
    t < zero(𝒯) ? one(𝒯) : zero(𝒯)
end

#--------------------------------------

"""Leaching model with exponential loss over time. The internal fields are

$(TYPEDFIELDS)

The constructor is

    function ExponentialLeaching(;
        λ::𝒯=1.0,
        C::𝒯=1.0,
        σ::𝒯=0.0,
        seed::Union{Nothing,Integer}=nothing,
    )

## Keyword Arguments

* The `λ` parameter is the exponential decay/loss rate
* `C` is a floor, above which the leache fraction will never cross.
* The `σ` parameter adds noise to the model, but should be used with caution because the noise is not symmetrical and can change the *average* leached fraction at a given point in time.
* If there is noise, the internal `rng` can be initialized using the `seed` keyword argument.
"""
struct ExponentialLeaching{𝒯} <: AbstractLeachingModel{𝒯}
    λ::𝒯
    C::𝒯
    σ::𝒯
    rng::Xoshiro
end

function ExponentialLeaching(;
    λ::𝒯=1.0,
    C::𝒯=1.0,
    σ::𝒯=0.0,
    seed::Union{Nothing,Integer}=nothing,
) where {𝒯}
    @assert λ >= zero(𝒯)
    @assert zero(𝒯) <= C <= one(𝒯)
    @assert σ >= zero(𝒯)
    ExponentialLeaching(λ, C, σ, Xoshiro(seed))
end

function (w::ExponentialLeaching)(t::𝒯)::𝒯 where {𝒯}
    if t < zero(𝒯)
        return one(𝒯)
    end
    @unpack λ, C, σ, rng = w
    𝓁 = C - C * exp(-λ * t)
    leachnoise(rng, 𝓁, σ)
end

#--------------------------------------

"""Leaching model with an average of exponential losses over time. The internal fields are

$(TYPEDFIELDS)

The constructor is

    function MultiExponentialLeaching(;
        λ::NTuple{𝒩,𝒯}=1.0,
        C::𝒯=1.0,
        σ::𝒯=0.0,
        seed::Union{Nothing,Integer}=nothing,
    )

## Keyword Arguments

* The `λ` parameter is tuple of decay/loss rates
* `C` is a floor, above which the leache fraction will never cross.
* The `σ` parameter adds noise to the model, but should be used with caution because the noise is not symmetrical and can change the *average* leached fraction at a given point in time.
* If there is noise, the internal `rng` can be initialized using the `seed` keyword argument.
"""
struct MultiExponentialLeaching{𝒯,𝒩} <: AbstractLeachingModel{𝒯}
    λ::NTuple{𝒩,𝒯}
    C::𝒯
    σ::𝒯
    rng::Xoshiro
end

function MultiExponentialLeaching(;
    λ::NTuple{𝒩,𝒯}=1.0,
    C::𝒯=1.0,
    σ::𝒯=0.0,
    seed::Union{Nothing,Integer}=nothing,
) where {𝒩,𝒯}
    @assert all(λ .>= 0)
    @assert 0 <= C <= 1
    @assert σ >= zero(𝒯)
    MultiExponentialLeaching(λ, C, σ, Xoshiro(seed))
end

function (w::MultiExponentialLeaching{𝒯,𝒩})(t::𝒯)::𝒯 where {𝒩,𝒯}
    if t < zero(𝒯)
        return one(𝒯)
    end
    @unpack λ, C, σ, rng = w
    𝓁 = C - C * sum(x -> exp(-x * t), λ) / 𝒩
    leachnoise(rng, 𝓁, σ)
end

#--------------------------------------
export integrand

"""Leaching model with seasonally varying loss over time, representing an annual cycle. The period of the cycles is exactly one year. The internal fields are

$(TYPEDFIELDS)

The constructor is

    function SeasonalLeaching(;
        λ=1.0,
        C=1.0,
        floor=0.0,
        power::Integer=1,
        phase=π,
        σ=0.0,
        seed::Union{Nothing,Integer}=nothing,
    )

## Keyword Arguments

* The `λ` parameter is the exponential decay/loss rate
* `C` is a floor, above which the leache fraction will never cross.
* `floor` sets the minimum value of the sinusoid that is integrated to produce monotonically increasing loss fractions. The details are not that important. It sets the minimum leaching *rate*. For example, if `floor` is zero, there is a point in each annual cycle where the leaching rate is also zero. As `floor` is increased, the minimum leaching rate over each annual cycle will increase.
* `power` defines how wide the cycles are, which is analogous to how long the "winter" is when the leaching rate is slow.
* `phase` shifts the cycle in time.
* The `σ` parameter adds noise to the model, but should be used with caution because the noise is not symmetrical and can change the *average* leached fraction at a given point in time.
* If there is noise, the internal `rng` can be initialized using the `seed` keyword argument.
"""
struct SeasonalLeaching{𝒯} <: AbstractLeachingModel{𝒯}
    λ::𝒯
    C::𝒯
    A::𝒯
    γ::UInt8
    ϕ::𝒯
    σ::𝒯
    rng::Xoshiro
end

function SeasonalLeaching(;
    λ=1.0,
    C=1.0,
    floor=0.0,
    power::Integer=1,
    phase=π,
    σ=0.0,
    seed::Union{Nothing,Integer}=nothing,
)
    @assert λ >= 0
    @assert 0 <= C < 1
    @assert 0 <= floor <= 1
    @assert power ∈ (1, 2, 3) "the power (exponent) can only be 1, 2, or 3"
    @assert σ >= 0
    λ, C, floor, phase, σ = promote(λ, C, floor, phase, σ)
    SeasonalLeaching(λ, C, 1 - floor, power |> UInt8, phase, σ, Xoshiro(seed))
end

function integrand(w::SeasonalLeaching, t)
    @unpack λ, C, A, γ, ϕ = w
    A * ((1 + cos(2π * t + ϕ)) / 2)^γ + (1 - A)
end

function (w::SeasonalLeaching{𝒯})(t::𝒯)::𝒯 where {𝒯}
    if t < zero(𝒯)
        return one(𝒯)
    end
    @unpack λ, C, A, γ, ϕ, σ, rng = w
    #! format: off
    integral = if γ == 1
        A * sin(π * t) * cos(π * t + ϕ) / (2π) - A * t / 2 + t
    elseif γ == 2
        (A * (8 * sin(2π * t + ϕ) + sin(2 * (2π * t + ϕ)) + 12π * t - 8 * sin(ϕ) - sin(2 * ϕ))) / (32π) - A * t + t
    elseif γ == 3
        (A * (45 * sin(2π * t + ϕ) + 9 * sin(2 * (2π * t + ϕ)) + sin(3 * (2π * t + ϕ)) + 60π * t - 45 * sin(ϕ) - 9 * sin(2 * ϕ) - sin(3 * ϕ))) / (192 * π) - A * t + t
    end
    #! format: on
    𝓁 = C - C * exp(-λ * integral)
    leachnoise(rng, 𝓁, σ)
end

end
