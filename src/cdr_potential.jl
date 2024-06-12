module CDRPotential

using PeriodicTable: elements
using Unitful: @u_str, ustrip

using DocStringExtensions

const VALID_MOBILE_CATIONS = [:Ca, :Mg, :Na, :K]

const GROUP = Dict{Symbol,Int}(:Ca => 2, :Mg => 2, :Na => 1, :K => 1)

const CO2_atomic_mass = 44.01u"u"

#------------------------------------------------------------------------------
export cationCO2, cdrpotential

function checkvalid(element::Symbol)::Nothing
    @assert element ∈ VALID_MOBILE_CATIONS "Element $element is not a valid mobile cation. Mobile cations are restricted to $VALID_MOBILE_CATIONS."
    nothing
end

"""$(TYPEDSIGNATURES)
Computes the mass of CO2 captured per unit mass of a mobile cation. The cation must be one of $VALID_MOBILE_CATIONS."""
function cationCO2(cation::Symbol)
    checkvalid(cation)
    GROUP[cation] * CO2_atomic_mass / elements[cation].atomic_mass
end

cationCO2(cation::String) = cation |> Symbol |> cationCO2

"""$(TYPEDSIGNATURES)
Computes the maximum amount of CDR available for a given feedstock composition."""
function cdrpotential(cation, concentration)
    cationCO2(cation) * concentration
end

"""$(TYPEDSIGNATURES)
Computes CDR potential from a named tuple of feedstock concentrations."""
function cdrpotential(conc::NamedTuple{𝒦,NTuple{𝒩,𝒯}})::𝒯 where {𝒦,𝒩,𝒯}
    sum(𝒦) do cation
        cdrpotential(cation, conc[cation])
    end
end

"""$(TYPEDSIGNATURES)
Computes CDR potential from a dict of feedstock concentrations."""
function cdrpotential(conc::Dict)
    sum(cdrpotential(k, conc[k]) for k ∈ keys(conc))
end

end
