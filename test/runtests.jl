using Pkg
using MultivariateDiscretization
using Test
using CSV
using DataFrames
using Distributions
using StatsBase
using Discretizers
using Statistics
using MultivariateStats
#push!(LOAD_PATH,"C:/Users/Niclas Popp/.julia/dev/MultivariateDiscretization.jl")


# number of dimensions
dims=30
# number of points
points=150


# generate test dataset
A = rand(Float64, (dims,dims))
Σ = A*A'
μ = rand(dims)
d1 = MvNormal(μ, 1/4*Σ)
Testdata = rand(d1,points)

bb = BayesianBlocks(Testdata,genes)
dr = DoaneRule(Testdata,genes)
cpd = CPD(Testdata, genes, cells, 3, :km, :manual, 20, :none)
ipd = greedy_IPD(Testdata,genes,35,:km)


@testset "MultivariateDiscretization.jl" begin

    @test typeof(ipd) == DataFrame
    @test typeof(cpd) == DataFrame
    @test typeof(bb) == DataFrame
    @test typeof(dr) == DataFrame

    cjs_ipd = CJS_empirical(Matrix(float.(ipd)),genes,cells,200)
    cjs_cpd = CJS_empirical(Matrix(float.(cpd)),genes,cells,200)
    cjs_bb = CJS_empirical(Matrix(float.(bb)),genes,cells,200)
    cjs_dr = CJS_empirical(Matrix(float.(dr)),genes,cells,200)
    cjs_or = CJS_empirical(Testdata,genes,cells,200)

    @test typeof(cjs_ipd) == Array{Float64,2}
    @test typeof(cjs_cpd) == Array{Float64,2}
    @test typeof(cjs_bb) == Array{Float64,2}
    @test typeof(cjs_dr) == Array{Float64,2}

end
