module Checks

using Distributions: UnivariateDistribution
using Unitful: Quantity, Units

#------------------------------------------------------------------------------
export ispositive, isnonnegative, isfractional, islessthanone

function ispositive(x::Quantity{𝒯}, u::Units)::Bool where {𝒯<:Real}
    (x > zero(𝒯) * u) ? true : false
end

function ispositive(x::𝒯, ::Units)::Bool where {𝒯<:Real}
    (x > zero(𝒯)) ? true : false
end

function ispositive(x::𝒯)::Bool where {𝒯<:Real}
    (x > zero(𝒯)) ? true : false
end

#--------------------------------------

function isnonnegative(x::Quantity{𝒯}, u::Units)::Bool where {𝒯<:Real}
    (x >= zero(𝒯) * u) ? true : false
end

function isnonnegative(x::𝒯, ::Units)::Bool where {𝒯<:Real}
    (x >= zero(𝒯)) ? true : false
end

function isnonnegative(x::𝒯)::Bool where {𝒯<:Real}
    (x >= zero(𝒯)) ? true : false
end

function isnonnegative(x::UnivariateDistribution)
    minimum(x) >= (x |> eltype |> zero) ? true : false
end

#--------------------------------------

function isfractional(x::Quantity{𝒯}, u::Units)::Bool where {𝒯<:Real}
    (zero(𝒯) * u <= x <= one(𝒯) * u) ? true : false
end

function isfractional(x::𝒯, ::Units)::Bool where {𝒯<:Real}
    (zero(𝒯) <= x <= one(𝒯)) ? true : false
end

function isfractional(x::𝒯)::Bool where {𝒯<:Real}
    (zero(𝒯) <= x <= one(𝒯)) ? true : false
end

#--------------------------------------

function islessthanone(x::Quantity{𝒯}, u::Units)::Bool where {𝒯<:Real}
    x < one(𝒯) * u
end

function islessthanone(x::𝒯, ::Units)::Bool where {𝒯<:Real}
    x < one(𝒯)
end

function islessthanone(x::𝒯)::Bool where {𝒯<:Real}
    x < one(𝒯)
end

end
