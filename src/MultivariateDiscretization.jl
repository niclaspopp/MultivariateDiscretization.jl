module MultivariateDiscretization

using DataFrames
using CSV
using StatsBase
using Discretizers
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
include("benchmark_measures.jl")

export
    greedy_IPD,
    CPD,
    CPD_clustered,
    BayesianBlocks,
    DoaneRule,
    CJS_empirical,
    UW,
    greedy_IPD_clustered

end
