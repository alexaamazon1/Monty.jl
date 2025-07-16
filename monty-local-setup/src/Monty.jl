module Monty

export ncore, nsample, nanalyte, analytes

function ncore end
function nsample end
function nanalyte end
function analytes end

include("checks.jl")
include("jitters.jl")
include("stencils.jl")
include("sampling.jl")
include("leaching_models.jl")
include("mixing_model.jl")
include("measuring.jl")
include("gaussian.jl")
include("simulating.jl")
include("cdr_potential.jl")

using Reexport: @reexport

@reexport using .Jitters
@reexport using .Stencils
@reexport using .Sampling
@reexport using .LeachingModels
@reexport using .MixingModel
@reexport using .Measuring
@reexport using .Gaussian
@reexport using .Simulating
@reexport using .CDRPotential

end
