
"""
    function modeConvergence(X, stops, numModes::Int)
Modal convergence check based on l2-norm of modes. The array stops contains
the ranges to investigate where stops[end] is used for reference. Subtraction of
mean is done in place.
"""
function modeConvergence(X, stops::AbstractArray{<: AbstractRange}, numModes::Int; subtractmean::Bool = false, method = :svd)

    if subtractmean
        X .-= mean(X,2)
    end

    numPODs = size(stops,1)

    # The POD to which all values are compared with
    if method == :svd
        maxPOD = PODsvd(X[:,stops[end]])[1].modes[:,1:numModes]
    elseif method == :eig
        maxPOD = PODeig(X[:,stops[end]])[1].modes[:,1:numModes]
    else
        error("supported methods, :eig, :svd")
    end


    output = zeros(Float64, numModes, numPODs)
    for i = 1:numPODs-1

        if method == :svd
            podRes = PODsvd(X[:,stops[i]])[1].modes[:,1:numModes]
        elseif method == :eig
            podRes = PODeig(X[:,stops[i]])[1].modes[:,1:numModes]
        end

        for j = 1:numModes

            # normalize the mode before comparing
            maxPODnorm = norm(maxPOD[:,j])
            podResnorm = norm(podRes[:,j])
            maxComp = maxPOD[:,j]/maxPODnorm
            podComp = podRes[:,j]/podResnorm

            # modes do not have a specific sign, compare the positive and negative
            # to find minimum.
            l2norm1 = norm(maxComp - podComp)
            l2norm2 = norm(maxComp + podComp)
            l2norm = minimum([l2norm1, l2norm2])
            output[j,i] = l2norm
        end

    end

    return output
end


"""
    function modeConvergence(X, W, stops, numModes::Int)
Same as function modeConvergence(X, stops, numModes::Int) but using weights.
"""
function modeConvergence(X, W::Vector, stops::AbstractArray{<: AbstractRange}, numModes::Int; subtractmean::Bool = false, method = :svd)

    if subtractmean
        X .-= mean(X,2)
    end

    numPODs = size(stops,1)

    # The POD to which all values are compared with
    if method == :svd
        maxPOD = PODsvd(X[:,stops[end]],W)[1].modes[:,1:numModes]
    elseif method == :eig
        maxPOD = PODeig(X[:,stops[end]],W)[1].modes[:,1:numModes]
    else
        error("supported methods, :eig, :svd")
    end


    output = zeros(Float64, numModes, numPODs)
    for i = 1:numPODs-1

        if method == :svd
            podRes = PODsvd(X[:,stops[i]],W)[1].modes[:,1:numModes]
        elseif method == :eig
            podRes = PODeig(X[:,stops[i]],W)[1].modes[:,1:numModes]
        end

        for j = 1:numModes

            # normalize the mode before comparing
            maxPODnorm = norm(maxPOD[:,j])
            podResnorm = norm(podRes[:,j])
            maxComp = maxPOD[:,j]/maxPODnorm
            podComp = podRes[:,j]/podResnorm

            # modes do not have a specific sign, compare the positive and negative
            # to find minimum.
            l2norm1 = norm(maxComp - podComp)
            l2norm2 = norm(maxComp + podComp)
            l2norm = minimum([l2norm1, l2norm2])
            output[j,i] = l2norm
        end

    end

    return output
end



"""
    function modeConvergence!(fileLoc, dataName, stops, numModes::Int)
Same as modeConvergence but data is loaded from a JLD file each iteration and
the POD is done in place.
"""
function modeConvergence!(fileLoc, dataName, stops::AbstractArray{<: AbstractRange}, numModes::Int; subtractmean::Bool = false, method = :svd)

    X = readdlm(fileLoc, ',')
    X = X[dataName]
    if subtractmean
        X .-= mean(X,2)
    end

    numPODs = size(stops,1)

    # The POD to which all values are compared with
    if method == :svd
        maxPOD = PODsvd(X[:,stops[end]])[1].modes[:,1:numModes]
    elseif method == :eig
        maxPOD = PODeig(X[:,stops[end]])[1].modes[:,1:numModes]
    else
        error("supported methods, :eig, :svd")
    end


    output = zeros(Float64, numModes, numPODs)
    for i = 1:numPODs-1
        X = readdlm(fileLoc, ',')
        X = X[dataName]
        if subtractmean
            X .-= mean(X,2)
        end

        if method == :svd
            podRes = PODsvd(X[:,stops[i]])[1].modes[:,1:numModes]
        elseif method == :eig
            podRes = PODeig(X[:,stops[i]])[1].modes[:,1:numModes]
        end

        for j = 1:numModes

            # normalize the mode before comparing
            maxPODnorm = norm(maxPOD[:,j])
            podResnorm = norm(podRes[:,j])
            maxComp = maxPOD[:,j]/maxPODnorm
            podComp = podRes[:,j]/podResnorm

            # modes do not have a specific sign, compare the positive and negative
            # to find minimum.
            l2norm1 = norm(maxComp - podComp)
            l2norm2 = norm(maxComp + podComp)
            l2norm = minimum([l2norm1, l2norm2])
            output[j,i] = l2norm
        end

    end

    return output
end



"""
    function modeConvergence!(fileLoc, dataName, W, stops, numModes::Int)
Same as modeConvergence!(fileLoc, dataName, stops, numModes::Int) but using weights.
"""
function modeConvergence!(fileLoc, dataName, W::Vector, stops::AbstractArray{<: AbstractRange}, numModes::Int; subtractmean::Bool = false, method = :svd)

    X = readdlm(fileLoc, ',')
    X = X[dataName]
    if subtractmean
        X .-= mean(X,2)
    end

    numPODs = size(stops,1)

    # The POD to which all values are compared with
    if method == :svd
        maxPOD = PODsvd(X[:,stops[end]],W)[1].modes[:,1:numModes]
    elseif method == :eig
        maxPOD = PODeig(X[:,stops[end]],W)[1].modes[:,1:numModes]
    else
        error("supported methods, :eig, :svd")
    end


    output = zeros(Float64, numModes, numPODs)
    for i = 1:numPODs-1
        X = readdlm(fileLoc, ',')
        X = X[dataName]
        if subtractmean
            X .-= mean(X,2)
        end

        if method == :svd
            podRes = PODsvd(X[:,stops[i]],W)[1].modes[:,1:numModes]
        elseif method == :eig
            podRes = PODeig(X[:,stops[i]],W)[1].modes[:,1:numModes]
        end

        for j = 1:numModes

            # normalize the mode before comparing
            maxPODnorm = norm(maxPOD[:,j])
            podResnorm = norm(podRes[:,j])
            maxComp = maxPOD[:,j]/maxPODnorm
            podComp = podRes[:,j]/podResnorm

            # modes do not have a specific sign, compare the positive and negative
            # to find minimum.
            l2norm1 = norm(maxComp - podComp)
            l2norm2 = norm(maxComp + podComp)
            l2norm = minimum([l2norm1, l2norm2])
            output[j,i] = l2norm
        end

    end

    return output
end
