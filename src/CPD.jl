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
# k: number of eigenvalues to keep
# d: discretizer used per principal components: :bb for bayesian blocks :uw for uniform width :km for kmeans
# bintype: type of binning, only used if disc=:uw, select :manual for manual setting of number of bins
# bins: number of bins, only used if bintype=:manual
function CPD(M, ndim::Int64, npoints::Int64, k::Int64, d=:km, bintype=:manual, bins=Int(round(ndim)),processing=:none)
    results = DataFrame()

    #calculate correlation matrix
    C=cor(M');

    #select principal components
    K = size(C, 1)
    ef=eigen(C)
    ef.values[K-k+1:K]
    S=ef.vectors[:, K-k+1:K]
    S=Matrix(S)

    #project onto principal components
    D=S'*M

    #determine cutpoints
    B = []
    for m in 1:k
        if d == :km
            R = kmeans(copy(convert(Array,D[m,:])'),bins)
            cutpoints=zeros(bins+1)
            cutpoints[1]=minimum(D[m,:])
            cutpoints[bins+1]=maximum(D[m,:])
            centers = sort(R.centers,dims=2)
            [cutpoints[i] = (centers[i-1]+centers[i])/2 for i in 2:bins]
            push!(B,cutpoints)
        else
            if d == :bb
                disc = LinearDiscretizer(binedges(DiscretizeBayesianBlocks(), D[m,:]))
            else
                if bintype==:manual
                    nbins=bins
                else
                    nbins=get_nbins(bintype, D[m,:])
                end

                disc = LinearDiscretizer(binedges(DiscretizeUniformWidth(nbins), D[m,:]))
            end
            cutpoints = zeros(nlabels(disc)+1)
            cutpoints[1] = extrema(disc,1)[1]
            #cutpoints = zeros(nlabels(disc))
            [cutpoints[i+1]=extrema(disc,i)[2] for i in 1:nlabels(disc)]
            push!(B,cutpoints)
        end
    end

    #reprojection the cutpoints
    cp_proj = []
    for i in B
        for j in i
            for l in 1:k
                push!(cp_proj, j*S[:,l])
            end
        end
    end

    projections = zeros(ndim,size(cp_proj)[1])

    for i in 1:ndim
        for j in 1:size(cp_proj)[1]
            projections[i,j]=cp_proj[j][i]
        end
        lindisc = LinearDiscretizer(sort(projections[i,:]))
        results[!,i] = [encode(lindisc,one) for one in M[i,:]]
    end

    if processing==:Integer
        results = postprocessing_integer(results,ndim,npoints)
    elseif processing==:Cutpoints
        results = postprocessing_cutpoints(results,ndim,npoints,projections)
    end

    return(results)
end

function CPD_clustered(M, ndim::Int64, npoints::Int64, c::Int64, p::Float64, d=:km)
    results = zeros(size(M))
    SPC = 1/2*(I-corspearman(M'))
    clust = hclust(SPC,linkage=:ward_presquared)
    clustered_data = cutree(clust;k=c)

    for i in 1:c
        selector = clustered_data .== i
        data = M[selector,:]
        k = Int(ceil(size(data)[1]/p))
        bins = get_nbins(:doane, data)
        res = CPD(data,size(data)[1],npoints,k,d,:manual,bins)
        #print("\n")
        #print(typeof(res))
        results[selector,:] .= Array(res)'
    end

    return(DataFrame(results))
end
