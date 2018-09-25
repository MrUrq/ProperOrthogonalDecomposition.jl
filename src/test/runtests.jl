using POD
using JLD
using Base.Test

# Define test matrix X with dimensions n×m, where n is number of data poitns and
# m is the number of snapshots
# W is the weight matrix representing the cell volume

X = rand(srand(1),1000,10)
X += 10
W = rand(srand(1),1:0.2:10,1000)

@testset "Standard POD" begin

    PODbase, Σ = POD(X)
    PODbaseEig, Σeig = PODeig(X)
    PODbaseSvd, Σsvd = PODsvd(X)

    meanPODbaseEig, meanΣeig = PODeig(copy(X), subtractmean = true)
    meanPODbaseSvd, meanΣsvd = PODsvd(copy(X), subtractmean = true)

    Σ₁ = 1050.362168606664
    Σ₂ = 8.301301177004774
    meanΣ₁ = 9.878480670625322
    meanΣ₂ = 8.30589586913458

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
        @test meanPODbaseEig.modes*meanPODbaseEig.coefficients ≈ X .- mean(X,2)
        @test meanPODbaseSvd.modes*meanPODbaseSvd.coefficients ≈ X .- mean(X,2)
    end
end

@testset "Weighted POD" begin

    PODbase, Σ = POD(X, W)
    PODbaseEig, Σeig = PODeig(X,W)
    PODbaseSvd, Σsvd = PODsvd(X,W)

    meanPODbaseEig, meanΣeig = PODeig(copy(X), W, subtractmean = true)
    meanPODbaseSvd, meanΣsvd = PODsvd(copy(X), W, subtractmean = true)

    Σ₁ = 2457.097163629827
    Σ₂ = 19.250259794561828
    meanΣ₁ = 23.33111514963416
    meanΣ₂ = 19.253445808771126

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
    A = [   0.010197212000627108 0.007070856617510225 0.004345509205119795 0
            0.6179627074745591 0.5476731024048103 0.12030902309667643 0
            1.377162205777348 0.8636384160114929 0.3197534098872419 0   ]
    Amean =[0.6322816398505274 0.5534565187905752 0.13120270384723295 0
            1.3821479511006756 0.8593667661723656 0.3083628970324381 0
            1.3536726312063883 1.1047589549024917 0.9957138077649103 0  ]
    W₂ = ones(1000)

    @testset "Convergence of number of included snapshots" begin
        @test_throws ErrorException modeConvergence(X,[1:4,1:6,1:8,1:10],3,method = :sdv)
        @test_throws ErrorException modeConvergence!(Pkg.dir("POD")*"/test/testdata.jld","X",[1:4,1:6,1:8,1:10],3,method = :sdv)

        @test modeConvergence(X,[1:4,1:6,1:8,1:10],3,method = :eig) ≈ A
        @test modeConvergence(X,[1:4,1:6,1:8,1:10],3,method = :svd) ≈ A
        @test modeConvergence!(Pkg.dir("POD")*"/test/testdata.jld","X",[1:4,1:6,1:8,1:10],3,method = :eig) ≈ A
        @test modeConvergence!(Pkg.dir("POD")*"/test/testdata.jld","X",[1:4,1:6,1:8,1:10],3,method = :svd) ≈ A

        @test modeConvergence(copy(X),[1:4,1:6,1:8,1:10],3,method = :eig,subtractmean = true) ≈ Amean
        @test modeConvergence(copy(X),[1:4,1:6,1:8,1:10],3,method = :svd,subtractmean = true) ≈ Amean
        @test modeConvergence!(Pkg.dir("POD")*"/test/testdata.jld","X",[1:4,1:6,1:8,1:10],3,method = :eig,subtractmean = true) ≈ Amean
        @test modeConvergence!(Pkg.dir("POD")*"/test/testdata.jld","X",[1:4,1:6,1:8,1:10],3,method = :svd,subtractmean = true) ≈ Amean
    end

    @testset "Convergence of number of included snapshots with one weights" begin
        @test_throws ErrorException modeConvergence(X,W₂,[1:4,1:6,1:8,1:10],3,method = :sdv)
        @test_throws ErrorException modeConvergence!(Pkg.dir("POD")*"/test/testdata.jld","X",[1:4,1:6,1:8,1:10],3,method = :sdv)

        @test modeConvergence(X,W₂,[1:4,1:6,1:8,1:10],3,method = :eig) ≈ A
        @test modeConvergence(X,W₂,[1:4,1:6,1:8,1:10],3,method = :svd) ≈ A
        @test modeConvergence!(Pkg.dir("POD")*"/test/testdata.jld","X",W₂,[1:4,1:6,1:8,1:10],3,method = :eig) ≈ A
        @test modeConvergence!(Pkg.dir("POD")*"/test/testdata.jld","X",W₂,[1:4,1:6,1:8,1:10],3,method = :svd) ≈ A

        @test modeConvergence(copy(X),W₂,[1:4,1:6,1:8,1:10],3,method = :eig,subtractmean = true) ≈ Amean
        @test modeConvergence(copy(X),W₂,[1:4,1:6,1:8,1:10],3,method = :svd,subtractmean = true) ≈ Amean
        @test modeConvergence!(Pkg.dir("POD")*"/test/testdata.jld","X",W₂,[1:4,1:6,1:8,1:10],3,method = :eig,subtractmean = true) ≈ Amean
        @test modeConvergence!(Pkg.dir("POD")*"/test/testdata.jld","X",W₂,[1:4,1:6,1:8,1:10],3,method = :svd,subtractmean = true) ≈ Amean
    end

    B = [   0.013795717122127843 0.008763206763625376 0
            1.076343280399526 0.5948830138034233 0
            1.349157219627426 1.2987842173266941 0  ]

    @testset "Convergence of step size" begin
        @test modeConvergence(X,[1:4:10,1:2:10,1:1:10],3) ≈ B
        @test modeConvergence!(Pkg.dir("POD")*"/test/testdata.jld","X",[1:4:10,1:2:10,1:1:10],3) ≈ B
    end
end
