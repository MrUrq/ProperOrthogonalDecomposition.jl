using ProperOrthogonalDecomposition
using Test
using StableRNGs
using Statistics
using DelimitedFiles

# Define test matrix X with dimensions n×m, where n is number of data poitns and
# m is the number of snapshots
# W is the weight matrix representing the cell volume

rng = StableRNGs.StableRNG(1)
X = rand(rng,1000,10)
X .+= 10
rng = StableRNGs.StableRNG(1)
W = rand(rng,1:0.2:10,1000)

@testset "Standard POD" begin

    PODbase, Σ = POD(X)
    PODbaseEig, Σeig = PODeigen(X)
    PODbaseSvd, Σsvd = PODsvd(X)

    meanPODbaseEig, meanΣeig = PODeigen(copy(X), subtractmean = true)
    meanPODbaseSvd, meanΣsvd = PODsvd(copy(X), subtractmean = true)

    Σ₁ = 1050.0370061587332
    Σ₂ = 8.413650551681984
    meanΣ₁ = 9.964071096830972
    meanΣ₂ = 8.42358075601378

    @testset "POD using eigenvalue decomposition" begin
        @test Σeig[1] ≈ Σ₁
        @test Σeig[end] ≈ Σ₂
        @test meanΣeig[1] ≈ meanΣ₁
        @test meanΣeig[end-1] ≈ meanΣ₂

    end

    @testset "POD using SVD" begin
        @test Σsvd[1] ≈ Σ₁
        @test Σsvd[end] ≈ Σ₂
        @test meanΣsvd[1] ≈ meanΣ₁
        @test meanΣsvd[end-1] ≈ meanΣ₂
    end

    @testset "Default method" begin
        @test maximum(abs.(abs.(PODbaseSvd.coefficients) .- abs.(PODbase.coefficients))) ≈ 0 atol=1e-8
        @test maximum(abs.(abs.(PODbaseSvd.modes) .- abs.(PODbase.modes))) ≈ 0 atol=1e-8
        @test maximum(abs.(Σsvd.-Σ)) ≈ 0 atol=1e-8
    end

    @testset "Method equality of SVD and Eig" begin
        @test maximum(abs.(abs.(PODbaseSvd.coefficients) .- abs.(PODbaseEig.coefficients))) ≈ 0 atol=1e-8
        @test maximum(abs.(abs.(PODbaseSvd.modes) .- abs.(PODbaseEig.modes))) ≈ 0 atol=1e-8
        @test maximum(abs.(abs.(meanPODbaseSvd.coefficients) .- abs.(meanPODbaseEig.coefficients))) ≈ 0 atol=1e-7
        @test maximum(abs.(abs.(meanPODbaseSvd.modes[:,1:end-1]) .- abs.(meanPODbaseEig.modes[:,1:end-1]))) ≈ 0 atol=1e-7
        @test maximum(abs.(Σeig.-Σsvd)) ≈ 0 atol=1e-8
    end

    @testset "Rebuild solution" begin
        @test PODbaseEig.modes*PODbaseEig.coefficients ≈ X
        @test PODbaseSvd.modes*PODbaseSvd.coefficients ≈ X
        @test meanPODbaseEig.modes*meanPODbaseEig.coefficients ≈ X .- mean(X,dims=2)
        @test meanPODbaseSvd.modes*meanPODbaseSvd.coefficients ≈ X .- mean(X,dims=2)
    end
end

@testset "Weighted POD" begin

    PODbase, Σ = POD(X, W)
    PODbaseEig, Σeig = PODeigen(X,W)
    PODbaseSvd, Σsvd = PODsvd(X,W)

    meanPODbaseEig, meanΣeig = PODeigen(copy(X), W, subtractmean = true)
    meanPODbaseSvd, meanΣsvd = PODsvd(copy(X), W, subtractmean = true)

    Σ₁ = 2462.447359829031
    Σ₂ = 19.570713247455668
    meanΣ₁ = 23.470360231176947
    meanΣ₂ = 19.595332053746322

    @testset "POD using eigenvalue decomposition" begin
        @test Σeig[1] ≈ Σ₁
        @test Σeig[end] ≈ Σ₂
        @test meanΣeig[1] ≈ meanΣ₁
        @test meanΣeig[end-1] ≈ meanΣ₂

    end

    @testset "POD using SVD" begin
        @test Σsvd[1] ≈ Σ₁
        @test Σsvd[end] ≈ Σ₂
        @test meanΣsvd[1] ≈ meanΣ₁
        @test meanΣsvd[end-1] ≈ meanΣ₂
    end

    @testset "Default method with weights" begin
        @test maximum(abs.(abs.(PODbaseSvd.coefficients) .- abs.(PODbase.coefficients))) ≈ 0 atol=1e-8
        @test maximum(abs.(abs.(PODbaseSvd.modes) .- abs.(PODbase.modes))) ≈ 0 atol=1e-8
        @test maximum(abs.(Σsvd.-Σ)) ≈ 0 atol=1e-8
    end

    @testset "Method equality of SVD and Eig with weights" begin
        @test maximum(abs.(abs.(PODbaseSvd.coefficients) .- abs.(PODbaseEig.coefficients))) ≈ 0 atol=1e-8
        @test maximum(abs.(abs.(PODbaseSvd.modes) .- abs.(PODbaseEig.modes))) ≈ 0 atol=1e-8
        @test maximum(abs.(Σeig.-Σsvd)) ≈ 0 atol=1e-8
    end

    @testset "Rebuild solution with weights" begin
        @test PODbaseEig.modes*PODbaseEig.coefficients ≈ X
        @test PODbaseSvd.modes*PODbaseSvd.coefficients ≈ X
    end
end

@testset "Mode convergence" begin
    A = [0.010164609073873544 0.006921747306424847 0.004321454495431182 0.0;
         0.7274414547632297 0.6102769856096316 0.2273767387120779 0.0;
         1.1710785171906266 1.1955870815537555 0.12812520121296264 0.0]
    W₂ = ones(1000)

    testdataPath = joinpath(dirname(pathof(ProperOrthogonalDecomposition)),"..","test","testdata.csv")

    @testset "Convergence of number of included snapshots" begin
        
        @test modeConvergence(X,PODeigen,[1:4,1:6,1:8,1:10],3) ≈ A
        @test modeConvergence(X,PODsvd,[1:4,1:6,1:8,1:10],3) ≈ A
        @test modeConvergence!(()->readdlm(testdataPath, ','),PODeigen!,[1:4,1:6,1:8,1:10],3) ≈ A
        @test modeConvergence!(()->readdlm(testdataPath, ','),PODsvd!,[1:4,1:6,1:8,1:10],3) ≈ A

    end

    @testset "Convergence of number of included snapshots with one weights" begin
        
        @test modeConvergence(X,x->PODeigen(x,W₂),[1:4,1:6,1:8,1:10],3) ≈ A
        @test modeConvergence(X,x->PODsvd(x,W₂),[1:4,1:6,1:8,1:10],3) ≈ A
        @test modeConvergence!(()->readdlm(testdataPath, ','),x->PODeigen!(x,W₂),[1:4,1:6,1:8,1:10],3) ≈ A
        @test modeConvergence!(()->readdlm(testdataPath, ','),x->PODsvd!(x,W₂),[1:4,1:6,1:8,1:10],3) ≈ A

    end

end
