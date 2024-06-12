module Measuring

using Random: AbstractRNG, Xoshiro, randn
using Unitful: @u_str
using UnPack: @unpack

using ..Checks
using ..MixingModel: Sample

import Base: show, getindex, zero

using DocStringExtensions

#------------------------------------------------------------------------------
export noise

"""$(TYPEDSIGNATURES)
Apply normally distributed noise to a number `x`, where the noise magnitude is expressed by the relative standard devation (`σᵣ`) of the number (the coefficient of variation or cv)."""
function noise(rng::AbstractRNG, x::𝒯, σᵣ::𝒯) where {𝒯<:Number}
    ϵ = x * (σᵣ * randn(rng, 𝒯))
    x + ϵ
end

#------------------------------------------------------------------------------
export Measurement, measure, zero

"""A `Measurement` behaves much like a [`Sample`](@ref) but has no core count. It represents observed values for a sample, after measurement noise is applied, and is produced by the [`measure`](@ref) methods.

The fields are
$(FIELDS)"""
struct Measurement{𝒩,𝒦,𝒯}
    concentrations::NamedTuple{𝒦,NTuple{𝒩,𝒯}}
    mass::𝒯
end

# indexing into the sample passes directly through to the concentrations
getindex(meas::Measurement, idx) = getindex(meas.concentrations, idx)

function zero(::Measurement{𝒩,𝒦,𝒯}) where {𝒩,𝒦,𝒯}
    Measurement(ntuple(i -> zero(𝒯), 𝒩) |> NamedTuple{𝒦}, zero(𝒯))
end

function show(io::IO, meas::Measurement{𝒩,𝒦,𝒯}) where {𝒩,𝒦,𝒯}
    @unpack mass = meas
    a = 𝒩 > 1 ? "analytes" : "analyte"
    println(io, "Measurement{$𝒯} ($𝒩 $a)\nmass = $mass")
    L = map(length ∘ string, 𝒦) |> maximum
    for (n, k) ∈ 𝒦 |> enumerate
        stick = (n < 𝒩) ? "├─" : "└─"
        println(io, "$stick $(rpad(k, L)) = $(meas[k])")
    end
end

#--------------------------------------

"""$(TYPEDSIGNATURES)
Computes measured analyte concentrations from a `Sample` by applying random normally distributed measurement noise to the exact values in the sample. Noise magnitude for analyte concentrations is defined by the relative standard deviation `σᵣ`, a named tuple with one value for each analyte in the sample. The measured sample mass noise is defined by relative standard deviation `σᵣ`. Returns a [`Measurement`](@ref) instance."""
function measure(
    rng::AbstractRNG,
    sample::Sample{𝒩,𝒦},
    σᵣ::NamedTuple{𝒦,NTuple{𝒩,𝒯}},
    σₘ::𝒯,
) where {𝒩,𝒦,𝒯<:Real}
    #apply measurement noise to concentrations
    conc = ntuple(𝒩) do j
        y = noise(rng, sample[j], σᵣ[j])
        @assert ispositive(y, u"kg/kg")
        return y
    end |> NamedTuple{𝒦}
    #include mass with it's own noise
    Measurement(conc, noise(rng, sample.mass, σₘ))
end

end
