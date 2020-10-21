using MultivariateDiscretization
using Test
# using DataFrames
# using Distributions
# using StatsBase
# using Statistics
# using MultivariateStats

#
# # number of dimensions
# dims=30
# # number of points
# points=150
#
#
# # generate test dataset
# A = rand(Float64, (dims,dims))
# Σ = A*A'
# μ = rand(dims)
# d1 = MvNormal(μ, 1/4*Σ)
# Testdata = rand(d1,points)
#
# bb = MultivariateDiscretization.BayesianBlocks(Testdata,dims)
# dr = MultivariateDiscretization.DoaneRule(Testdata,dims)
# cpd = MultivariateDiscretization.CPD(Testdata, dims, points, 3, :km, :manual, 20, :none)
# ipd = MultivariateDiscretization.greedy_IPD(Testdata,dims,35,:km)
#
#
# @testset "MultivariateDiscretization.jl" begin
#
#     @test typeof(ipd) == DataFrame
#     @test typeof(cpd) == DataFrame
#     @test typeof(bb) == DataFrame
#     @test typeof(dr) == DataFrame
#
#     cjs_ipd = CJS_empirical(Matrix(float.(ipd)),dims,points,200)
#     cjs_cpd = CJS_empirical(Matrix(float.(cpd)),dims,points,200)
#     cjs_bb = CJS_empirical(Matrix(float.(bb)),dims,points,200)
#     cjs_dr = CJS_empirical(Matrix(float.(dr)),dims,points,200)
#     cjs_or = CJS_empirical(Testdata,dims,points,200)
#
#     @test typeof(cjs_ipd) == Array{Float64,2}
#     @test typeof(cjs_cpd) == Array{Float64,2}
#     @test typeof(cjs_bb) == Array{Float64,2}
#     @test typeof(cjs_dr) == Array{Float64,2}
#
# end
