
"""
    PODsvd(X)

Uses the SVD method of snapshots to calculate the POD basis of X. Method of
snapshots is efficient when number of data points n > number of snapshots m.
"""
function PODsvd(X; subtractmean::Bool = false)

    Xcop = deepcopy(X)
    PODsvd!(Xcop, subtractmean = subtractmean)

end

"""
    PODsvd(X,W)

Same as `PODsvd(X)` but uses weights for each data point. The weights are equal to the
cell volume for a volume mesh.
"""
function PODsvd(X,W::AbstractVector; subtractmean::Bool = false)

    Xcop = deepcopy(X)
    PODsvd!(Xcop, W, subtractmean = subtractmean)

end

"""
    PODsvd!(X)

Same as `PODsvd(X)` but overwrites memory.
"""
function PODsvd!(X; subtractmean::Bool = false)

    if subtractmean
        X .-= mean(X,2)
    end

    # Economy sized SVD
    F = svdfact!(X)

    # Mode coefficients
    a = Diagonal(F.S)*F.Vt

    POD = PODBasis(a, F.U)

    return POD, F.S

end

"""
    PODsvd!(X,W)

Same as `PODsvd!(X)` but uses weights for each data point. The weights are equal to the
cell volume for a volume mesh.
"""
function PODsvd!(X,W::AbstractVector; subtractmean::Bool = false)

    if subtractmean
        X .-= mean(X,2)
    end

    # Take into account the weights
    X = sqrt.(W).*X

    # Economy sized SVD
    F = svdfact!(X)

    # Mode coefficients
    a = Diagonal(F.S)*F.Vt
    phi = 1 ./sqrt.(W).*F.U

    POD = PODBasis(a, phi)

    return POD, F.S

end
