function BayesianBlocks(M,ndim)

    results = DataFrame()

    for i in 1:ndim

        if length(unique(M[i,:]))<3
            res_temp = zeros(length(M[i,:]))

            for u in 1:length(unique(M[i,:]))
                res_temp[M[i,:] .== unique(M[i,:])[u]] .= u
            end

            results[!,i] = res_temp
        else

            disc = LinearDiscretizer(binedges(DiscretizeBayesianBlocks(), M[i,:]))
            results[!,i] = [encode(disc,one) for one in M[i,:]]

        end

    end

    return(results)

end
