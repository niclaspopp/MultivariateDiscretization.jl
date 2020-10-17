module MultivariateDiscretization

using DataFrames
using CSV
using StatsBase
using Discretizers
using Plots
using Distributions
using Statistics
using LinearAlgebra
using MultivariateStats
using FreqTables
using Clustering
using InformationMeasures
using Distances

include("DoaneRule.jl")
include("BayesianBlocks.jl")
include("CPD.jl")
include("IPD.jl")
include("UniformWidth.jl")
include("postprocessing.jl")
include("benchmark_measures.jl")

export
    greedy_IPD,
    greedy_IPD_cutpoints,
    CPD,
    BayesianBlocks,
    DoaneRule,
    CJS_empirical,
    UW,
    greedy_IPD_clustered

end
