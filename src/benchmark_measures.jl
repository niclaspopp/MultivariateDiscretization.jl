function CJS_empirical(F::Matrix,ndim::Int64,nsamples::Int64,nsteps::Int64)

    results = zeros(ndim,ndim)

    for k in 1:ndim
        for l in k:ndim
            lbound = minimum(F[:,k])
            rbound = maximum(F[:,k])

            #print(edf(F[:,k],rbound))

            int=0.0

            for u in 1:nsteps

                x = lbound+(u-1)*(rbound-lbound)/nsteps
                p=edf(F[:,k],x)
                q=edf(F[:,l],x)

                int+=(p*log2(2*p/(p+q))+1/2 * 1/log(2) * (q-p))/nsteps

            end

            results[k,l] = int

        end
    end

    return(results+results')
end

function edf(data::Array{Float64},x::Float64)
    return(length(data[data .<= x])/length(data))
end
