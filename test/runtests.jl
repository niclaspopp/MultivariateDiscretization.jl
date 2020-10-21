using MultivariateDiscretization
using Test
using DataFrames
using Distributions
using StatsBase
using Statistics
using LinearAlgebra
using MultivariateStats


# number of dimensions
dims=10
# number of points
points=100


# generate test dataset
A = rand(Float64, (dims,dims))
Σ = A*A'
μ = rand(dims)
#d1 = MvNormal(μ, 1/4*Σ)
d1 = MvNormal(μ, I)
Testdata = rand(d1,points)


@testset "MultivariateDiscretization.jl" begin

    bb = MultivariateDiscretization.BayesianBlocks(Testdata,dims)
    dr = MultivariateDiscretization.DoaneRule(Testdata,dims)
    cpd = MultivariateDiscretization.CPD_clustered(Testdata,dims,points,3,float(3))
    ipd = MultivariateDiscretization.greedy_IPD_clustered(Testdata,dims,3,3)
    uw = MultivariateDiscretization.UW(Testdata,dims,3)

    @test typeof(ipd) == DataFrame
    @test typeof(cpd) == DataFrame
    @test typeof(bb) == DataFrame
    @test typeof(dr) == DataFrame
    @test typeof(uw) == DataFrame

    cjs_ipd = MultivariateDiscretization.CJS_empirical(Matrix(float.(ipd)),dims,points,200)
    cjs_cpd = MultivariateDiscretization.CJS_empirical(Matrix(float.(cpd)),dims,points,200)
    cjs_bb = MultivariateDiscretization.CJS_empirical(Matrix(float.(bb)),dims,points,200)
    cjs_dr = MultivariateDiscretization.CJS_empirical(Matrix(float.(dr)),dims,points,200)
    cjs_uw = MultivariateDiscretization.CJS_empirical(Matrix(float.(uw)),dims,points,200)
    cjs_or = MultivariateDiscretization.CJS_empirical(Testdata,dims,points,200)

    @test typeof(cjs_ipd) == Array{Float64,2}
    @test typeof(cjs_cpd) == Array{Float64,2}
    @test typeof(cjs_bb) == Array{Float64,2}
    @test typeof(cjs_dr) == Array{Float64,2}
    @test typeof(cjs_uw) == Array{Float64,2}
    @test typeof(cjs_or) == Array{Float64,2}

end
