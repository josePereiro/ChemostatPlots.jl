module ChemostatPlots

import Chemostat
const Ch = Chemostat
const ChU = Chemostat.Utils
import LinearAlgebra: normalize
using Plots
using Distributions

include("plot_marginal.jl")
# include("plot_xi.jl")
# include("plot_stoi_err_beta.jl")
# include("plot_stoi_err_xi.jl")
# include("plot_mean_stoi_err_beta.jl")
# include("plot_mean_stoi_err_xi.jl")
# include("plot_norm_stoi_err_beta.jl")

end
