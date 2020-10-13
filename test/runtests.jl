push!(LOAD_PATH, pwd())
using MultivariateDiscretization
using Test
using CSV
using DataFrames
using Distributions
using StatsBase
using Discretizers
using Statistics
using MultivariateStats


# number of genes for checking
genes=30
# number of cells (first column contains gene names)
cells=150


####  CELLDATASET   ####
celldata_complete= CSV.read("test/data/norm_counts_hvg.csv", DataFrame)
gene_names = celldata_complete[1:genes,1]
celldata=celldata_complete[1:genes,1:cells+1]
norm_type = :log
Mp = Matrix(celldata[:,2:size(celldata,2)])
if norm_type == :log
    Mp=map(x->log.(1+x),Mp)
end
dt = fit(UnitRangeTransform, Mp, dims=1, unit=true)
M = StatsBase.transform(dt, Mp)
Testdata = Mp .- mean(Mp, dims = 2)

ipd = greedy_IPD(Testdata,genes,20,:km)
cpd = CPD(Testdata,genes,cells,3,:km,:manual, 5, :Integer)
bb = BayesianBlocks(Testdata,genes)
dr = DoaneRule(Testdata,genes)


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
