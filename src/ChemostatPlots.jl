module ChemostatPlots

import Chemostat
const Ch = Chemostat
const ChU = Chemostat.Utils
const ChLP = Chemostat.LP
import LinearAlgebra: normalize
using Plots
using Distributions

include("plot_marginal.jl")
include("plot_projection2D.jl")

end
