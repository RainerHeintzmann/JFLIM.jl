module JFLIM

using Statistics, BioformatsLoader, CUDA, NDTools, FourierTools, IndexFunArrays
using InverseModeling
using Plots 
using MicroscopyTools # for soft_delta
# using DeconvOptim # for conv_aux
using ChainRulesCore
# using LinearAlgebra # for norm

# from the DeconvOptim toolbox:
include("conv.jl")
include("forward_models.jl")

export loss_gaussian, loss_anscombe_pos, loss_anscombe, loss_poisson

include("flim_load.jl")
include("flim_fit.jl")
include("show_fit.jl")

end # module JFLIM