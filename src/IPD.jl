function interaction_distance(M::Array{Float64},data1::DataFrame,data2::DataFrame,ndim::Int64)

    k=length(data1[1])
    l=length(data2[1])

    sum = 0.0
    maxs = zeros(ndim)
    for i in 1:ndim
        maxs[i] = maximum(M[i,:])
    end

    for one in 1:k
        for two in 1:k
            prod=1.0
            for i in 1:ndim
                prod=prod*(maxs[i]-max(data1[one,i],data1[two,i]))
            end
            sum=sum+(1/k^2)*prod

        end
    end

    for one in 1:k
        for two in 1:l
            prod=1.0
            for i in 1:ndim
                prod=prod*(maxs[i]-max(data1[one,i],data2[two,i]))
            end
            sum=sum-(2/(k*l))*prod
        end
    end

    for one in 1:l
        for two in 1:l
            prod=1.0
            for i in 1:ndim
                prod=prod*(maxs[i]-max(data2[one,1],data2[two,1]))
            end
            sum=sum+(1/l^2)*prod
        end
    end

    return(abs(sum)^(1/2))
end

# dim_i: current dimension
function get_IDS_largeI(M::Array{Float64},B::DataFrame,labels::Array{Int64},ndim::Int64,dim_i::Int64)

    ids = zeros(length(labels)-1)
    M_temp = deepcopy(M[1:end .!= dim_i, 1:end])

    for i in 1:length(labels)-1

        data1_arr=[]
        data2_arr=[]

        for j in 1:ndim-1
            push!(data1_arr,M_temp[j,:][B[:,1] .== labels[i]])
            push!(data2_arr,M_temp[j,:][B[:,1] .== labels[i+1]])
        end

        data1 = DataFrame(data1_arr)
        data2 = DataFrame(data2_arr)

        ids[i]=interaction_distance(M_temp,data1,data2,ndim-1)

    end

    return(ids)

end

function largeI(bin::Array{Int64},interaction_distances::Array{Float64}, threshold::Float64)

    sum = 0

    if length(bin)>1
        for x in 1:length(bin)-1
            if interaction_distances[bin[x]]>=threshold
                sum = sum+1
            end
        end
    end

    return(sum)

end

function determine_threshold(interaction_distances::Array{Float64},T::Int64)

    tert=max(1,Int(round((T-1)/20)))

    return(sort(interaction_distances)[tert])

end

function logstar(x::Int64)

    y=deepcopy(x)
    res = 0.0

    while log2(y)>0
        y=log2(y)
        res = res+y
    end

    return(res)
end

function LN(z::Int64)
    return(logstar(z)+log2(2.8654))
end

function encode_disc(T::Int64,k::Int64)
    return(LN(k)+log2(binomial(T-1,k-1)))
end


function Lbid(merges::Array{Any},T::Int64)

    sum = 0
    for bin in merges

        sum = sum+LN(length(bin))-(1+length(bin))*log2(length(bin)/T)

    end

    return(sum)
end

function Lmh(merges::Array{Any},interaction_distances::Array{Float64},threshold::Float64)

    sum = 0.0
    I = 0

    for bin in merges

        I = largeI(bin,interaction_distances,threshold)
        if I>0
            sum = sum + LN(I) + I*log2(length(bin)-1)
        end

    end

    return(sum)

end

function encode_discrete_data(merges::Array{Any},interaction_distances::Array{Float64},T::Int64,threshold::Float64)
    return(Lbid(merges,T)+Lmh(merges,interaction_distances,threshold))
end

function encode_error(merges::Array{Any})
    sum=0
    for bin in merges
        sum = sum+length(bin)*log2(length(bin))
    end

    return(sum)

end


# M: original data
# B: Micro binned data
# merges: array of array that form the macro bins
# k: number of macro bins
# t: threshold
# ndim: number of dimensions
# T: nnumber of micro bins
function practical_score(M::Array{Float64},B::DataFrame,merges::Array{Any}, k::Int64, threshold::Float64, ndim::Int64,T::Int64, interaction_distances::Array{Float64})

    return(encode_error(merges)+encode_discrete_data(merges,interaction_distances,T,threshold)+encode_disc(T,k))

end



#### --------------------------------------------------------------------------------------------------------------------------------------------------####
#### --------------------------------------------------------------------------------------------------------------------------------------------------####

# Input for greedy_IPD_v2:
#M: data
# nm: number of micro bins
# ndim: number of dimension
function greedy_IPD(M::Array{Float64},ndim::Int64,T::Int64,disc=:km)

    results = DataFrame()
    #cuts = DataFrame()
    #results_arr=[]

    for i in 1:ndim

        if mod(i,100)==0
            print("\n")
            print("dimension ")
            print(i)
        end

        micro_binned = DataFrame();

        # label for results in dataframe
        #label = string("dimension"," ",string(i))

        # start by uniform discretization for micro bins, seperate per dimension !!!

        if length(unique(M[i,:]))<=5

            res_temp = zeros(length(M[i,:]))

            for u in 1:length(unique(M[i,:]))
                res_temp[M[i,:] .== unique(M[i,:])[u]] .= u
            end

            results[!,i] = res_temp

        else

            if length(unique(M[i,:]))<=T
                T=length(unique(M[i,:]))
            end

            if disc==:km
                R = kmeans(vec(M[i,:])',T)
                cutpoints=zeros(T+1)
                cutpoints[1]=minimum(M[i,:])
                cutpoints[T+1]=maximum(M[i,:])
                centers = sort(R.centers,dims=2)
                [cutpoints[i] = (centers[i-1]+centers[i])/2 for i in 2:T]

                # discretize the data
                lindisc = LinearDiscretizer(sort(cutpoints))
                micro_binned[!, :1] = [encode(lindisc,one) for one in M[i,:]]
            else

                lindisc = LinearDiscretizer(binedges(DiscretizeUniformWidth(T), M[i,:]))
                micro_binned[!, :1] = [encode(lindisc,one) for one in M[i,:]]

            end


            #get the lables
            labels_micro = sort(unique(micro_binned[:,:1]))

            #initialize number of macro bins
            k=T


            # determine the interaction distances
            int_dists = get_IDS_largeI(M,micro_binned,labels_micro,ndim,i)

            # determine threshlod for largeI
            t=determine_threshold(int_dists,T)

            merges=[]
            for j in labels_micro
                push!(merges,[j])
            end


            scores=[]
            # practical score at the beginning
            MDL_ref=practical_score(M,micro_binned,merges,k-1,t,ndim,T,int_dists)
            append!(scores,MDL_ref)

            while minimum(scores)<=MDL_ref

                optimal_merge = 0

                if length(merges)>1
                    for j in 1:length(merges)-1

                        temp_merges = deepcopy(merges)
                        append!(temp_merges[j],temp_merges[j+1])
                        deleteat!(temp_merges, j+1)

                        mdl = practical_score(M,micro_binned,temp_merges,k-1,t,ndim,T,int_dists)

                        # check if the current tuple is the best so far
                        if mdl < minimum(scores)
                            optimal_merge = j
                            push!(scores,mdl)
                        end

                    end
                end

                if optimal_merge != 0

                    append!(merges[optimal_merge],merges[optimal_merge+1])
                    deleteat!(merges, optimal_merge+1)
                    MDL_ref = minimum(scores)
                    k=k-1

                else
                    MDL_ref=minimum(scores)-1
                end

            end


            results_dim = zeros(length(M[i,:]))
            for index in 1:length(merges)
                macrobin = merges[index]
                for id in macrobin
                    results_dim[micro_binned[!,:1] .== id] .= index
                end

            end

            results[!,i] = results_dim
        end

    end

    return(results)

end




function greedy_IPD_Cutpoints(M::Array{Float64},ndim::Int64,T::Int64,disc=:km)

    results = DataFrame()
    cuts = DataFrame()
    #results_arr=[]

    for i in 1:ndim

        if mod(i,100)==0
            print("\n")
            print("dimension ")
            print(i)
        end

        micro_binned = DataFrame();

        # label for results in dataframe
        #label = string("dimension"," ",string(i))

        # start by uniform discretization for micro bins, seperate per dimension !!!

        if length(unique(M[i,:]))<=5

            res_temp = zeros(length(M[i,:]))

            for u in 1:length(unique(M[i,:]))
                res_temp[M[i,:] .== unique(M[i,:])[u]] .= unique(M[i,:])[u]
            end

            results[!,i] = res_temp

        else

            if length(unique(M[i,:]))<=T
                T=length(unique(M[i,:]))
            end

            if disc==:km
                R = kmeans(vec(M[i,:])',T)
                cutpoints=zeros(T+1)
                cutpoints[1]=minimum(M[i,:])
                cutpoints[T+1]=maximum(M[i,:])
                centers = sort(R.centers,dims=2)
                [cutpoints[i] = (centers[i-1]+centers[i])/2 for i in 2:T]

                # discretize the data
                lindisc = LinearDiscretizer(sort(cutpoints))
                micro_binned[!, :1] = [encode(lindisc,one) for one in M[i,:]]
            else

                lindisc = LinearDiscretizer(binedges(DiscretizeUniformWidth(T), M[i,:]))
                micro_binned[!, :1] = [encode(lindisc,one) for one in M[i,:]]

            end

            results[!,i] = cutpoints


            #get the lables
            labels_micro = sort(unique(micro_binned[:,:1]))

            #initialize number of macro bins
            k=T


            # determine the interaction distances
            int_dists = get_IDS_largeI(M,micro_binned,labels_micro,ndim,i)

            # determine threshlod for largeI
            t=determine_threshold(int_dists,T)

            merges=[]
            for j in labels_micro
                push!(merges,[j])
            end


            scores=[]
            # practical score at the beginning
            MDL_ref=practical_score(M,micro_binned,merges,k-1,t,ndim,T,int_dists)
            append!(scores,MDL_ref)

            while minimum(scores)<=MDL_ref

                optimal_merge = 0

                if length(merges)>1
                    for j in 1:length(merges)-1

                        temp_merges = deepcopy(merges)
                        append!(temp_merges[j],temp_merges[j+1])
                        deleteat!(temp_merges, j+1)

                        mdl = practical_score(M,micro_binned,temp_merges,k-1,t,ndim,T,int_dists)

                        # check if the current tuple is the best so far
                        if mdl < minimum(scores)
                            optimal_merge = j
                            push!(scores,mdl)
                        end

                    end
                end

                if optimal_merge != 0

                    append!(merges[optimal_merge],merges[optimal_merge+1])
                    deleteat!(merges, optimal_merge+1)
                    MDL_ref = minimum(scores)
                    k=k-1

                else
                    MDL_ref=minimum(scores)-1
                end

            end


            results_dim = zeros(length(M[i,:]))
            for index in 1:length(merges)
                macrobin = merges[index]
                maxpoint = maximum(macrobin)
                minpoint = minimum(macrobin)
                for id in macrobin
                    results_dim[micro_binned[!,:1] .== id] .= (cutpoints[maxpoint]-cutpoints[minpoint])/2
                end

            end

            results[!,i] = results_dim
        end

    end

    return(results)

end
