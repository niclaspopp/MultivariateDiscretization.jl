function postprocessing_integer(df::DataFrame,ndim::Int64,npoints::Int64)
    #res=DataFrame(undef, size(df))
    for i in 1:ndim
        ints = unique(df[:,i])
        for j in 1:length(ints)
            df[!,i][df[:,i] .== ints[j]] .= j
        end
    end
    return df
end

function postprocessing_cutpoints(df::DataFrame,ndim::Int64,npoints::Int64, cutpoints::Array{Float64})
    #res=DataFrame(undef, size(df))
    df=float.(df)
    for i in 1:ndim
        ints = unique(df[:,i])
        for j in 1:length(ints)
            df[!,i][df[:,i] .== ints[j]] .= .5*(sort(cutpoints[i,:])[j]-sort(cutpoints[i,:])[j+1])
        end
    end
    return df
end
