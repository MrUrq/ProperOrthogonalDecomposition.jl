var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#ProperOrthogonalDecomposition-1",
    "page": "Home",
    "title": "ProperOrthogonalDecomposition",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Introduction-1",
    "page": "Home",
    "title": "Introduction",
    "category": "section",
    "text": "ProperOrthogonalDecomposition is a Julia package for performing the Proper Orthogonal modal Decomposition (POD) technique. The technique has been used  to, among other things, extract turbulent flow features. The POD methods available in this package is the Singular Value Decomposition (SVD) based method and the eigen-decomposition based method of snapshots. The method is snapshots is the most commonly used method for fluid flow analysis where the number of  datapoints is larger than the number of snapshots.The POD technique goes under several names; Karhunen-Loèven (KL), Principal Component Analysis (PCA) and Hotelling analysis. The method has been used for error analysis, reduced order modeling, fluid flow reconstruction, turbulent flow feature extraction, etc. A descriptive overview of the method is given in [1]."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "The package is registered and can be installed with Pkg.add.julia> Pkg.add(\"ProperOrthogonalDecomposition\")"
},

{
    "location": "index.html#Reference-1",
    "page": "Home",
    "title": "Reference",
    "category": "section",
    "text": "[1]: Taira et al. \"Modal Analysis of Fluid Flows: An Overview\", arXiv:1702.01453 [physics], () http://arxiv.org/abs/1702.01453"
},

{
    "location": "man/POD.html#",
    "page": "POD",
    "title": "POD",
    "category": "page",
    "text": ""
},

{
    "location": "man/POD.html#POD-1",
    "page": "POD",
    "title": "POD",
    "category": "section",
    "text": "Each method returns a tuple containing the pod basis of type PODBasis{T} and  the corresponding singular values. The singular values are related to each modes importance to the dataset. "
},

{
    "location": "man/POD.html#ProperOrthogonalDecomposition.PODeigen-Tuple{Any}",
    "page": "POD",
    "title": "ProperOrthogonalDecomposition.PODeigen",
    "category": "method",
    "text": "PODeigen(X)\n\nUses the eigenvalue method of snapshots to calculate the POD basis of X. Method of snapshots is efficient when number of data points n > number of snapshots m.\n\n\n\n\n\n"
},

{
    "location": "man/POD.html#ProperOrthogonalDecomposition.PODeigen!-Tuple{Any}",
    "page": "POD",
    "title": "ProperOrthogonalDecomposition.PODeigen!",
    "category": "method",
    "text": "PODeigen!(X)\n\nSame as PODeigen(X) but overwrites memory.\n\n\n\n\n\n"
},

{
    "location": "man/POD.html#Method-of-snapshots-1",
    "page": "POD",
    "title": "Method of snapshots",
    "category": "section",
    "text": "The eigen-decomposition based method of snapshots is the most commonly used  method for fluid flow analysis where the number of datapoints is larger than the number of snapshots.PODeigen(X; subtractmean::Bool = false)PODeigen!(X; subtractmean::Bool = false)"
},

{
    "location": "man/POD.html#ProperOrthogonalDecomposition.PODsvd-Tuple{Any}",
    "page": "POD",
    "title": "ProperOrthogonalDecomposition.PODsvd",
    "category": "method",
    "text": "PODsvd(X)\n\nUses the SVD based decomposition technique to calculate the POD basis of X. \n\n\n\n\n\n"
},

{
    "location": "man/POD.html#ProperOrthogonalDecomposition.PODsvd!-Tuple{Any}",
    "page": "POD",
    "title": "ProperOrthogonalDecomposition.PODsvd!",
    "category": "method",
    "text": "PODsvd!(X)\n\nSame as PODsvd(X) but overwrites memory.\n\n\n\n\n\n"
},

{
    "location": "man/POD.html#Singular-Value-Decomposition-based-method-1",
    "page": "POD",
    "title": "Singular Value Decomposition based method",
    "category": "section",
    "text": "The SVD based approach is also available and is more robust against roundoff errors. PODsvd(X; subtractmean::Bool = false)PODsvd!(X; subtractmean::Bool = false)"
},

{
    "location": "man/POD.html#Example-1",
    "page": "POD",
    "title": "Example",
    "category": "section",
    "text": "Here we will artifically create data which is PODed and then extract the first mode.t, x = range(0, stop=30, length=50), range(-10, stop=30, length=120)\n\nXgrid = [i for i in x, j in t]\ntgrid = [j for i in x, j in t]\n\nf1 = sech.(Xgrid.-3.5) .* 10.0 .* cos.(0.5 .*tgrid)\nf2 = cos.(Xgrid) .* 1.0 .* cos.(2.5 .*tgrid)\nf3 = sech.(Xgrid.+5.0) .* 4.0 .* cos.(1.0 .*tgrid)\n\nY = f1+f2+f3\n\nusing PlotlyJS # hide\n\nfunction plotpoddata(Y) # hide\n    trace = surface(x=Xgrid,y=tgrid,z=Y,colorscale=\"Viridis\", cmax=7.5, cmin=-7.5) # hide\n    layout = Layout(height=440, # hide\n                    scene = (   xaxis=attr(title=\"Space\"), # hide\n                                yaxis=attr(title=\"Time\"), # hide\n                                zaxis=attr(title=\"z\",range=[-10,10])), # hide\n                    margin=attr(l=30, r=30, b=20, t=90), # hide\n                    ) # hide\n    plot(trace, layout) # hide\nend # hide\np = plotpoddata(Y) # hide\npkgpath = abspath(joinpath(dirname(Base.find_package(\"ProperOrthogonalDecomposition\")), \"..\")) # hide\nsavedir = joinpath(pkgpath,\"docs\",\"src\",\"assets\",\"poddata.html\") # hide\nPlotlyJS.savehtml(p,savedir,:embed) # hideOur data Y looks like this    <iframe src=\"../assets/poddata.html\" height=\"540\" width=\"765\" frameborder=\"0\" seamless=\"seamless\" scrolling=\"no\"></iframe>Now we POD the data and reconstruct the dataset using only the first mode.using ProperOrthogonalDecomposition # hide\nres, singularvals  = POD(Y)\nreconstructFirstMode = res.modes[:,1:1]*res.coefficients[1:1,:]\n\np = plotpoddata(reconstructFirstMode) # hide\npkgpath = abspath(joinpath(dirname(Base.find_package(\"ProperOrthogonalDecomposition\")), \"..\")) # hide\nsavedir = joinpath(pkgpath,\"docs\",\"src\",\"assets\",\"podfirstmode.html\") # hide\nPlotlyJS.savehtml(p,savedir,:embed) # hideNote that the above used POD(Y) which defaults to the SVD based apparoch. The first mode over the time series looks like this    <iframe src=\"../assets/podfirstmode.html\" height=\"540\" width=\"765\" frameborder=\"0\" seamless=\"seamless\" scrolling=\"no\"></iframe>"
},

{
    "location": "man/weightedPOD.html#",
    "page": "Weighted POD",
    "title": "Weighted POD",
    "category": "page",
    "text": ""
},

{
    "location": "man/weightedPOD.html#ProperOrthogonalDecomposition.PODeigen-Tuple{Any,AbstractArray{T,1} where T}",
    "page": "Weighted POD",
    "title": "ProperOrthogonalDecomposition.PODeigen",
    "category": "method",
    "text": "PODeigen(X,W)\n\nSame as PODeigen(X) but uses weights for each data point. The weights are equal to the cell volume for a volume mesh.\n\n\n\n\n\n"
},

{
    "location": "man/weightedPOD.html#ProperOrthogonalDecomposition.PODeigen!-Tuple{Any,AbstractArray{T,1} where T}",
    "page": "Weighted POD",
    "title": "ProperOrthogonalDecomposition.PODeigen!",
    "category": "method",
    "text": "PODeigen!(X,W)\n\nSame as PODeigen!(X) but uses weights for each data point. The weights are equal to the cell volume for a volume mesh.\n\n\n\n\n\n"
},

{
    "location": "man/weightedPOD.html#ProperOrthogonalDecomposition.PODsvd-Tuple{Any,AbstractArray{T,1} where T}",
    "page": "Weighted POD",
    "title": "ProperOrthogonalDecomposition.PODsvd",
    "category": "method",
    "text": "PODsvd(X,W)\n\nSame as PODsvd(X) but uses weights for each data point. The weights are equal to the cell volume for a volume mesh.\n\n\n\n\n\n"
},

{
    "location": "man/weightedPOD.html#ProperOrthogonalDecomposition.PODsvd!-Tuple{Any,AbstractArray{T,1} where T}",
    "page": "Weighted POD",
    "title": "ProperOrthogonalDecomposition.PODsvd!",
    "category": "method",
    "text": "PODsvd!(X,W)\n\nSame as PODsvd!(X) but uses weights for each data point. The weights are equal to the cell volume for a volume mesh.\n\n\n\n\n\n"
},

{
    "location": "man/weightedPOD.html#Weighted-POD-1",
    "page": "Weighted POD",
    "title": "Weighted POD",
    "category": "section",
    "text": "When performing the POD method it is assumed that the datapoints are equidistantly spaced.  This assumption makes the method sensistive to the local mesh resolution. To make the method mesh independent, a vector with weights for each datapoint can be supplied. Typically the weights are chosen to be the cell volume, although the face area can be used in the case of a plane. PODeigen(X,W::AbstractVector; subtractmean::Bool = false)PODeigen!(X,W::AbstractVector; subtractmean::Bool = false)PODsvd(X,W::AbstractVector; subtractmean::Bool = false)PODsvd!(X,W::AbstractVector; subtractmean::Bool = false)"
},

{
    "location": "man/weightedPOD.html#Example-1",
    "page": "Weighted POD",
    "title": "Example",
    "category": "section",
    "text": "Here we create the same data as in the previous example; however, we refine the  mesh locally, at x>7.5 && x<=30 and plot the reconstructed data from the first mode."
},

{
    "location": "man/weightedPOD.html#Non-uniform-grid-*without*-weights-1",
    "page": "Weighted POD",
    "title": "Non-uniform grid without weights",
    "category": "section",
    "text": "\nt, xcoarse = range(0, stop=30, length=50), range(-10, stop=7.5, length=30)\nxfine = range(7.5+step(xcoarse), stop=30, length=1000)\n\nx = [xcoarse...,xfine...]\nXgrid = [i for i in x, j in t]\ntgrid = [j for i in x, j in t]\n\nf1 = sech.(Xgrid.-3.5) .* 10.0 .* cos.(0.5 .*tgrid)\nf2 = cos.(Xgrid) .* 1.0 .* cos.(2.5 .*tgrid)\nf3 = sech.(Xgrid.+5.0) .* 4.0 .* cos.(1.0 .*tgrid)\n\nY = f1+f2+f3\n\nusing PlotlyJS # hide\nfunction plotpoddata(Y) # hide\n    trace = surface(x=Xgrid,y=tgrid,z=Y,colorscale=\"Viridis\", cmax=7.5, cmin=-7.5) # hide\n    layout = Layout(height=440, # hide\n                    scene = (   xaxis=attr(title=\"Space\"), # hide\n                                yaxis=attr(title=\"Time\"), # hide\n                                zaxis=attr(title=\"z\",range=[-10,10])), # hide\n                    margin=attr(l=30, r=30, b=20, t=90), # hide\n                    ) # hide\n    plot(trace, layout) # hide\nend # hide\np = plotpoddata(Y) # hide\n\nusing ProperOrthogonalDecomposition # hide\nres, singularvals  = POD(Y)\nreconstructFirstMode = res.modes[:,1:1]*res.coefficients[1:1,:]\np = plotpoddata(reconstructFirstMode) # hide\npkgpath = abspath(joinpath(dirname(Base.find_package(\"ProperOrthogonalDecomposition\")), \"..\")) # hide\nsavedir = joinpath(pkgpath,\"docs\",\"src\",\"assets\",\"finemeshfirstmode.html\") # hide\nPlotlyJS.savehtml(p,savedir,:embed) # hideAnd the first three singular values.singularvals[1:3]The first mode has changed due to the local mesh refinement compated to the previously presented case with equidistant mesh.    <iframe src=\"../assets/finemeshfirstmode.html\" height=\"540\" width=\"765\" frameborder=\"0\" seamless=\"seamless\" scrolling=\"no\"></iframe>"
},

{
    "location": "man/weightedPOD.html#Non-uniform-grid-*with*-weights-1",
    "page": "Weighted POD",
    "title": "Non-uniform grid with weights",
    "category": "section",
    "text": "Using the volume weighted formulation removes the mesh depedency and we get the correct modes back. grid_resolution = [repeat([step(xcoarse)],length(xcoarse));\n                   repeat([step(xfine)],length(xfine))]\nres, singularvals  = POD(Y,grid_resolution)\nreconstructFirstMode = res.modes[:,1:1]*res.coefficients[1:1,:]\n\np = plotpoddata(reconstructFirstMode) # hide\npkgpath = abspath(joinpath(dirname(Base.find_package(\"ProperOrthogonalDecomposition\")), \"..\")) # hide\nsavedir = joinpath(pkgpath,\"docs\",\"src\",\"assets\",\"finemeshfirstmodeweighted.html\") # hide\nPlotlyJS.savehtml(p,savedir,:embed) # hideAnd the first three singular values.singularvals[1:3]    <iframe src=\"../assets/finemeshfirstmodeweighted.html\" height=\"540\" width=\"765\" frameborder=\"0\" seamless=\"seamless\" scrolling=\"no\"></iframe>"
},

{
    "location": "man/weightedPOD.html#Uniform-grid-with-weights-1",
    "page": "Weighted POD",
    "title": "Uniform grid with weights",
    "category": "section",
    "text": "Compare the singular values from the above two cases with the singular values  from the weighted POD on the equidistant mesh.t, x = range(0, stop=30, length=50), range(-10, stop=30, length=120)\ngrid_resolution = repeat([step(x)],length(x))\n\nXgrid = [i for i in x, j in t]\ntgrid = [j for i in x, j in t]\n\nf1 = sech.(Xgrid.-3.5) .* 10.0 .* cos.(0.5 .*tgrid)\nf2 = cos.(Xgrid) .* 1.0 .* cos.(2.5 .*tgrid)\nf3 = sech.(Xgrid.+5.0) .* 4.0 .* cos.(1.0 .*tgrid)\n\nY = f1+f2+f3\n\nres, singularvals  = POD(Y,grid_resolution)\nnothing # hideAnd the first three singular values.singularvals[1:3]"
},

{
    "location": "man/convergence.html#",
    "page": "Mode convergence",
    "title": "Mode convergence",
    "category": "page",
    "text": ""
},

{
    "location": "man/convergence.html#ProperOrthogonalDecomposition.modeConvergence-Tuple{AbstractArray,Any,AbstractArray{#s11,N} where N where #s11<:AbstractRange,Int64}",
    "page": "Mode convergence",
    "title": "ProperOrthogonalDecomposition.modeConvergence",
    "category": "method",
    "text": "function modeConvergence(X::AbstractArray, PODfun, stops::AbstractArray{<: AbstractRange}, numModes::Int)\n\nModal convergence check based on l2-norm of modes. The array stops contains the ranges to investigate where stops[end] is used as the reference modes. The numModes largest modes are compared to reduce the computational time. The function used to POD the data is supplied through PODfun\n\n\n\n\n\n"
},

{
    "location": "man/convergence.html#ProperOrthogonalDecomposition.modeConvergence!-Tuple{Any,Any,AbstractArray{#s11,N} where N where #s11<:AbstractRange,Int64}",
    "page": "Mode convergence",
    "title": "ProperOrthogonalDecomposition.modeConvergence!",
    "category": "method",
    "text": "function modeConvergence!(loadFun, PODfun, stops::AbstractArray{<: AbstractRange}, numModes::Int)\n\nSame as modeConvergence(X::AbstractArray, PODfun, stops::AbstractArray{<: AbstractRange}, numModes::Int)  but here the data is reloaded for each comparision so that an inplace POD method can be used to reduce maximum memory usage.\n\n\n\n\n\n"
},

{
    "location": "man/convergence.html#Mode-convergence-1",
    "page": "Mode convergence",
    "title": "Mode convergence",
    "category": "section",
    "text": "Functionality to investigate convergence is supplied in this packages where the convergence in time and frequency can be investigated.modeConvergence(X::AbstractArray, PODfun, stops::AbstractArray{<: AbstractRange}, numModes::Int)modeConvergence!(loadFun, PODfun, stops::AbstractArray{<: AbstractRange}, numModes::Int)"
},

{
    "location": "man/convergence.html#Example-1",
    "page": "Mode convergence",
    "title": "Example",
    "category": "section",
    "text": ""
},

{
    "location": "man/convergence.html#Convergence-in-time-1",
    "page": "Mode convergence",
    "title": "Convergence in time",
    "category": "section",
    "text": "using PlotlyJS # hide\nusing ProperOrthogonalDecomposition\n\nt, x = range(0, stop=100, length=100), range(-10, stop=30, length=120)\n\nXgrid = [i for i in x, j in t]\ntgrid = [j for i in x, j in t]\n\nf1 = sech.(Xgrid.-3.5) .* 10.0 .* cos.(0.5 .*tgrid)\nf2 = cos.(Xgrid) .* 1.0 .* cos.(2.5 .*tgrid)\nf3 = sech.(Xgrid.+5.0) .* 4.0 .* cos.(1.0 .*tgrid)\n\nY = f1+f2+f3\n\n#Array of ranges we\'re interested in investigating\nranges = Array{UnitRange{Int64}}(undef,40)         \n\n#Ranges of interest starting from 3 timesteps \nsubset = range(3, stop=size(Y,2), length=length(ranges))\nfor i = 1:length(ranges)\n    ranges[i] = 1:round(Int,subset[i])\nend\n\nconvergence = modeConvergence(Y,PODeigen,ranges,3)\n\n\nfunction plotconvergence(subset,convergence) # hide\n    x=round.(Int,subset) # hide\n\n    trace1 = scatter(;x=x, y=convergence[1,:], # hide\n                        mode=\"markers\", name=\"Mode 1\", # hide\n                        marker_size=12) # hide\n\n    trace2 = scatter(;x=x, y=convergence[2,:], # hide\n                        mode=\"markers\", name=\"Mode 2\", # hide\n                        marker_size=12) # hide\n\n    trace3 = scatter(;x=x, y=convergence[3,:], # hide\n                        mode=\"markers\", name=\"Mode 3\", # hide\n                        marker_size=12) # hide\n    \n    data = [trace1, trace2, trace3] # hide\n    layout = Layout(height=440, # hide\n                    title=\"Time Convergence\", # hide\n                    xaxis=attr(title=\"Time\"), # hide\n                    yaxis=attr(title=\"Norm difference \"), # hide\n                    margin=attr(l=100, r=30, b=50, t=90), # hide\n                                ) # hide\n    plot(data, layout) # hide\nend # hide\n\np = plotconvergence(subset,convergence) # hide\npkgpath = abspath(joinpath(dirname(Base.find_package(\"ProperOrthogonalDecomposition\")), \"..\")) # hide\nsavedir = joinpath(pkgpath,\"docs\",\"src\",\"assets\",\"convergenceTime.html\") # hide\nPlotlyJS.savehtml(p,savedir,:embed) # hideThe history of convergence indicates the point at which additional data no longer provides additional  information to the POD modes.    <iframe src=\"../assets/convergenceTime.html\" height=\"540\" width=\"765\" frameborder=\"0\" seamless=\"seamless\" scrolling=\"no\"></iframe>"
},

{
    "location": "man/convergence.html#Convergence-inplace-1",
    "page": "Mode convergence",
    "title": "Convergence inplace",
    "category": "section",
    "text": "Datasets can quickly become large which is why an inplace method is available where the user supplies a function to load the data.using DelimitedFiles\n\n#Anonymous function with zero arguments\nloadFun = ()->readdlm(\"path/to/data/dataset.csv\", \',\')\n\n#POD the data inplace and reload it into memory each time.\nconvergence = modeConvergence!(loadFun,PODeigen!,ranges,3)This can also be done for a weighted POD withconvergence = modeConvergence!(loadFun,X->PODeigen!(X,W),ranges,3)note: Note\nThe use of a delimited files, such as a *.csv in the above example,  is not advisable if memory is a concern. Use a binary file format such as HDF5 for example. "
},

{
    "location": "man/convergence.html#Convergence-in-frequency-1",
    "page": "Mode convergence",
    "title": "Convergence in frequency",
    "category": "section",
    "text": "Just as we can investigate the time history needed for the mode to be converged,  we can also investigate the sampling frequency needed. This is done by supplying the  ranges as subsampled sets of the full time history.using PlotlyJS # hide\nusing ProperOrthogonalDecomposition\n\nt, x = range(0, stop=50, length=1000), range(-10, stop=30, length=120)\n\nXgrid = [i for i in x, j in t]\ntgrid = [j for i in x, j in t]\n\nf1 = sech.(Xgrid.-3.5) .* 10.0 .* cos.(0.5 .*tgrid)\nf2 = cos.(Xgrid) .* 1.0 .* cos.(2.5 .*tgrid)\nf3 = sech.(Xgrid.+5.0) .* 4.0 .* cos.(1.0 .*tgrid)\n\nY = f1+f2+f3\n\n#Array of ranges we\'re interested in investigating \nsubset = 100:-3:1 #Sub-sampling starts at every 100:th timestep\nranges = Array{StepRange{Int64,Int64}}(undef,length(subset))         \n\nfor i = 1:length(ranges)\n    ranges[i] = 1:round(Int,subset[i]):length(t)\nend\n\nconvergence = modeConvergence(Y,PODeigen,ranges,3)\n\n\nfunction plotconvergence(subset,convergence) # hide\n    x=1 ./((length(t)/last(t)) ./round.(Int,subset)) # hide\n\n    trace1 = scatter(;x=x, y=convergence[1,:], # hide\n                        mode=\"markers\", name=\"Mode 1\", # hide\n                        marker_size=12) # hide\n\n    trace2 = scatter(;x=x, y=convergence[2,:], # hide\n                        mode=\"markers\", name=\"Mode 2\", # hide\n                        marker_size=12) # hide\n\n    trace3 = scatter(;x=x, y=convergence[3,:], # hide\n                        mode=\"markers\", name=\"Mode 3\", # hide\n                        marker_size=12) # hide\n    \n    data = [trace1, trace2, trace3] # hide\n    layout = Layout(height=440, # hide\n                    title=\"Sampling Frequency Convergence\", # hide\n                    xaxis=attr(title=\"1/Freq.\"), # hide\n                    yaxis=attr(title=\"Norm difference \"), # hide\n                    margin=attr(l=100, r=30, b=50, t=90), # hide\n                                ) # hide\n    plot(data, layout) # hide\nend # hide\n\np = plotconvergence(subset,convergence) # hide\npkgpath = abspath(joinpath(dirname(Base.find_package(\"ProperOrthogonalDecomposition\")), \"..\")) # hide\nsavedir = joinpath(pkgpath,\"docs\",\"src\",\"assets\",\"convergenceFreq.html\") # hide\nPlotlyJS.savehtml(p,savedir,:embed) # hidenote: Note\nThe data point where 1/f = 1.25 indicates that Mode 2 and Mode 3 are far from converged, this sudden jump is likely due to the relative importance of the modes switching at this sampling frequency. This does not necessarily mean that the  modes themselves are poorly represented.    <iframe src=\"../assets/convergenceFreq.html\" height=\"540\" width=\"765\" frameborder=\"0\" seamless=\"seamless\" scrolling=\"no\"></iframe>"
},

]}
