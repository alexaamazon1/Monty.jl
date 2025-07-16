module Stencils

using Random: Xoshiro
using Meshes: Point, Vec
using UnPack: @unpack

using ..Jitters

import ..ncore

using DocStringExtensions

#------------------------------------------------------------------------------
export AbstractStencil

abstract type AbstractStencil{𝒯} end

"""$(TYPEDSIGNATURES)
Returns the number of cores in a stencil"""
ncore(S::AbstractStencil) = S.N

(S::AbstractStencil)(T::Tuple{Point{2},Int}) = S(T[1], T[2])

#--------------------------------------
export SingleCoreStencil

"""An empty type representing a single core. The SingleCoreStencil constructor has no arguments

    SingleCoreStencil()
"""
struct SingleCoreStencil{𝒯} <: AbstractStencil{𝒯} end

SingleCoreStencil(𝒯::Type) = SingleCoreStencil{𝒯}()

ncore(::SingleCoreStencil) = one(UInt16)

(S::SingleCoreStencil)(point::Point{2})::Point{2} = point

function (S::SingleCoreStencil)(point::Point{2}, core::Integer)::Point{2}
    @assert core == 1 "SingleCoreStencil expects core number equal to 1, but got $core"
    S(point)
end

show(io::IO, ::SingleCoreStencil) =
    println(io, "SingleCoreStencil\n└─ cores = 1")

#--------------------------------------
export RandomStencil

"""A core stencil representing randomly drawn core locations

    RandomStencil(N::Integer, args...; kwargs...)

Creates a `RandomStencil` with `N` cores, where the `args` and `kwargs` are passed through to the [`Jitter`](@ref) constructor"""
struct RandomStencil{𝒯,𝒟} <: AbstractStencil{𝒯}
    N::UInt16
    J::Jitter{𝒯,𝒟}
end

"""$(TYPEDSIGNATURES)
"""
function RandomStencil(N::Integer, args...; kwargs...)
    RandomStencil(N |> UInt16, Jitter(args...; kwargs...))
end

(S::RandomStencil)(point::Point{2})::Point{2} = point |> S.J

(S::RandomStencil)(point::Point{2}, ::Integer)::Point{2} = S(point)

function show(io::IO, S::RandomStencil{𝒯,𝒟}) where {𝒯,𝒟}
    println(
        io,
        "RandomStencil{$𝒯}\n├─ cores = $(S.N)\n└─ distribution = $(S.J.d)",
    )
end

#--------------------------------------
export CircleStencil

"""A core stencil representing a circle of cores

A CircleStencil is defined by the number of cores and the radius of the circle

    CircleStencil(N::Integer, r)"""
struct CircleStencil{𝒯<:Number} <: AbstractStencil{𝒯}
    N::UInt16
    r::𝒯
    s::Vector{𝒯}
    c::Vector{𝒯}
end

function CircleStencil(N, r::𝒯) where {𝒯}
    θ = ntuple(i -> convert(𝒯, 2π * (i - 1) / N), N)
    CircleStencil(N |> UInt16, r, sin.(θ) |> collect, cos.(θ) |> collect)
end

function (S::CircleStencil{𝒯})(
    point::Point{2},
    core::Integer,
)::Point{2} where {𝒯}
    @unpack N, r, s, c = S
    @assert 0 < core <= N "expected core number between 0 and $N, but got $core"
    point + Vec(r * c[core], r * s[core])
end

function show(io::IO, S::CircleStencil{𝒯}) where {𝒯}
    println(io, "CircleStencil{$𝒯}\n├─ cores = $(S.N)\n└─ radius = $(S.r)")
end

#--------------------------------------
export HubSpokeStencil

"""A core stencil representing a circle of cores with one point in the center of the circle.

    HubSpokeStencil(N::Integer, r)

The constructor takes the number of cores and a radius for the circle"""
struct HubSpokeStencil{𝒯} <: AbstractStencil{𝒯}
    ○::CircleStencil{𝒯}
end

HubSpokeStencil(N::Integer, r::Real) = HubSpokeStencil(CircleStencil(N - 1, r))

function (S::HubSpokeStencil)(point::Point{2}, core::Integer)::Point{2}
    (core == 1) ? point : S.○(point, core - 1)
end

ncore(S::HubSpokeStencil) = ncore(S.○) + 1

function show(io::IO, S::HubSpokeStencil{𝒯}) where {𝒯}
    println(
        io,
        "HubSpokeStencil{$𝒯}\n├─ cores = $(ncore(S))\n└─ radius = $(S.○.r)",
    )
end

#--------------------------------------
export LineStencil

"""A core stencil representing a line of cores.

    LineStencil(N::Integer, δ::𝒯, θ::𝒯) where {𝒯}

The constructor takes the number of cores, the distance between each core in the line (`δ`), and the angle of the line (`θ`)
"""
struct LineStencil{𝒯<:Number} <: AbstractStencil{𝒯}
    N::UInt16
    s::𝒯
    c::𝒯
    Δ::𝒯
end

function LineStencil(N::Integer, δ::𝒯, θ::𝒯) where {𝒯}
    s, c = sincos(θ)
    Δ = (N - 1) * δ
    LineStencil(N |> UInt16, s, c, Δ)
end

LineStencil(N::Integer, δ::𝒯) where {𝒯<:Real} = LineStencil(N, δ, zero(𝒯))

function (S::LineStencil{𝒯})(
    point::Point{2,𝒯},
    core::Integer,
)::Point{2} where {𝒯}
    @unpack N, s, c, Δ = S
    @assert 0 < core <= N
    f = convert(𝒯, (core - 1) / (N - 1))
    h = one(𝒯) / 2
    point + Vec(c * Δ * (f - h), s * Δ * (f - h))
end

function show(io::IO, S::LineStencil{𝒯}) where {𝒯}
    δ = S.Δ / (S.N - 1)
    θ = acos(S.c)
    println(
        io,
        "LineStencil{$𝒯}\n├─ cores = $(S.N)\n├─ spacing = $δ\n└─ angle = $θ",
    )
end

end
