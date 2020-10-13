function get_nbins(alg, data)

    n = length(data)

    if alg == :sqrt
        nbins = ceil(Int, sqrt(n))
    elseif alg == :sturges
        nbins = ceil(Int, log(2,n)) + 1
    elseif alg == :rice
        nbins = ceil(Int, 2*cbrt(n))
    elseif alg == :doane
        g₁ = moment(data, 3)
        σ = sqrt((6*(n-2))/((n+1)*(n+3)))
        nbins = ceil(Int, 1 + log(2,n) + log(2, 1+abs(g₁)/σ))
    elseif alg == :scott
        σ = std(data)
        binwidth = 3.5*σ/cbrt(n)
        lo, hi = extrema(data)
        nbins = ceil(Int, (hi - lo)/binwidth)
    elseif alg == :fd
        binwidth = 2*iqr(data)/cbrt(n)
        lo, hi = extrema(data)
        nbins = ceil((hi - lo)/binwidth)
    else
        binwidth = 2*iqr(data)/cbrt(n)
        lo, hi = extrema(data)
        nbins_fd = ceil(Int, (hi - lo)/binwidth)
        nbins_sturges = ceil(Int, log(2,n)) + 1
        nbins = max(nbins_fd, nbins_sturges)
    end

    return(round(nbins))
end

# M: data
# ndim: number of dimensions
# bintype: type of binning, select :manual for manual setting of number of bins
# bins: number of bins, only used if bintype=:manual
function DoaneRule(M,ndim)
    results = DataFrame()

    for i in 1:ndim

        if length(unique(M[i,:]))<5
            res_temp = zeros(length(M[i,:]))

            for u in 1:length(unique(M[i,:]))
                res_temp[M[i,:] .== unique(M[i,:])[u]] .= u
            end

            results[!,i] = res_temp
        else

            nbins = get_nbins(:doane, M[i,:])
            disc = LinearDiscretizer(binedges(DiscretizeUniformWidth(nbins), M[i,:]))
            results[!,i] = [encode(disc,one) for one in M[i,:]]

        end

    end

    return(results)

end
