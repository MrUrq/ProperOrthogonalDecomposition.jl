
"""
    PODeig(X)

Uses the eigenvalue method of snapshots to calculate the POD basis of X. Method of
snapshots is efficient when number of data points n > number of snapshots m.
"""
function PODeig(X; subtractmean::Bool = false)

    Xcop = deepcopy(X)
    PODeig!(Xcop, subtractmean = subtractmean)

end


"""
    PODeig(X,W)

Same as `PODeig(X)` but uses weights for each data point. The weights are equal to the
cell volume for a volume mesh.
"""
function PODeig(X,W::AbstractVector; subtractmean::Bool = false)

    Xcop = deepcopy(X)
    PODeig!(Xcop,W, subtractmean = subtractmean)

end

"""
    PODeig!(X,W)

Same as `PODeig!(X)` but uses weights for each data point. The weights are equal to the
cell volume for a volume mesh.
"""
function PODeig!(X,W::AbstractVector; subtractmean::Bool = false)

    if subtractmean
        X .-= mean(X,2)
    end

    # Number of snapshots
    m = size(X,2)

    # Correlation matrix for method of snapshots
    C = X'*Diagonal(W)*X

    # Eigen Decomposition
    E = eigfact!(C)
    eigVals = E.values
    eigVects = E.vectors

    # Sort the eigen vectors
    sortInd = sortperm(abs.(eigVals)/m,rev=true)
    eigVects = eigVects[:,sortInd]
    eigVals = eigVals[sortInd]

    # Diagonal matrix containing the square roots of the eigenvalues
    S = sqrt.(abs.(eigVals))

    # Construct the modes and coefficients
    phi = X*eigVects*Diagonal(1 ./S)
    a = Diagonal(S)*eigVects'

    POD = PODBasis(a, phi)

    return POD, S

end


"""
    PODeig!(X)

Same as `PODeig(X)` but overwrites memory.
"""
function PODeig!(X; subtractmean::Bool = false)

    if subtractmean
        X .-= mean(X,2)
    end

    # Number of snapshots
    m = size(X,2)

    # Correlation matrix for method of snapshots
    C = X'*X

    # Eigen Decomposition
    E = eigfact!(C)
    eigVals = E.values
    eigVects = E.vectors

    # Sort the eigen vectors
    sortInd = sortperm(abs.(eigVals)/m,rev=true)
    eigVects = eigVects[:,sortInd]
    eigVals = eigVals[sortInd]

    # Diagonal matrix containing the square roots of the eigenvalues
    S = sqrt.(abs.(eigVals))

    # Construct the modes and coefficients
    phi = X*eigVects*Diagonal(1 ./S)
    a = Diagonal(S)*eigVects'

    POD = PODBasis(a, phi)

    return POD, S

end
