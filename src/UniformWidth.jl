function UW(M,ndim,nbins)
    results = zeros(size(M))

    for i in 1:ndim

        if length(unique(M[i,:]))<3
            res_temp = zeros(length(M[i,:]))

            for u in 1:length(unique(M[i,:]))
                res_temp[M[i,:] .== unique(M[i,:])[u]] .= u
            end

            results[i,:] = res_temp
        else

            disc = LinearDiscretizer(binedges(DiscretizeUniformWidth(nbins), M[i,:]))
            results[i,:] = [encode(disc,one) for one in M[i,:]]

        end

    end

    return(DataFrame(results))

    return(results)

end
