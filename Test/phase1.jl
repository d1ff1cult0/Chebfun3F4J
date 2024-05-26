using LinearAlgebra
using FastChebInterp
using InteractiveUtils
using Statistics

# Include helper functions
folder_path = joinpath(@__DIR__, "Private")

files = readdir(folder_path)
for file in files
    include(joinpath(folder_path, file))
end

# Include data structures and related functions
include("chebfun3f.jl")

@time function phase1(f :: Function, enableLog::Bool = false, enableTime:: Bool = false, grid::Int64 = 17)
    # Informatie uit chebfunpref() van matlab
    pseudoLevel = eps() :: Float64

    if enableLog == 1
        println("\nChebfun3F constructor:")
        evalsold = 0;
    end

    #Initialize
    vectorize = false :: Bool
    n = fill(grid,3) :: Vector{Int64}
    r = fill(6,3) :: Vector{Int64}
    tol = pseudoLevel :: Float64
    reltol = 0.0 :: Float64
    # Setup cf3F
    cf3F = chebfun3F([[]], [[]], [[]], [[]], [], 0, 0) :: chebfun3F
    cheb(i,n) = 0 .- (cos.((i.-1).*(pi/(n-1))))
    reffun(n) = (floor.(Int64,sqrt(2).^(floor.(2*log2.(n)) .+ 1)) .+ 1)

    # Check fast vectorization
    try
        A = f(1:2,1:2,1:2) :: Vector{Float64}
        vectorize = true
        if isScalar(A)
            g(x,y,z) = f(x,y,z) + 0*x + 0*y + 0*z # ensure vectorization
            f = g
        end
    catch
        vectorize = false
    end

    t1 = 0
    t2 = 0
    t3 = 0

    #Main loop
    happy::Bool = false
    while ~happy
        # extra declarations
        I = Vector{Int64}() ; J = Vector{Int64}() ; K = Vector{Int64}()
        IT2 = Vector{Int64}() ; IT3 = Vector{Int64}()
        JT1 = Vector{Int64}() ; JT3 = Vector{Int64}()
        KT1 = Vector{Int64}() ; KT2 = Vector{Int64}()
        I2 = Vector{Int64}() ; J2 = Vector{Int64}() ; K2 = Vector{Int64}()
        Uc = Matrix{Float64}(undef,0,0) ; Vc = Matrix{Float64}(undef,0,0) ; Wc = Matrix{Float64}(undef,0,0)
        Uf = Matrix{Float64}(undef,0,0) ; Vf = Matrix{Float64}(undef,0,0) ; Wf = Matrix{Float64}(undef,0,0)
        
        # Phase 1
        happyPhase1 = false :: Bool
        while !happyPhase1
            J = initializeIndexRandomly(r[2], n[2])
            K = initializeIndexRandomly(r[3], n[3])
            # Handle to evaluate tensor entries of T_c
            n1 :: Int64, n2 :: Int64, n3 :: Int64 = n
            ff(i,j,k) = f(cheb(i,n1), cheb(j,n2), cheb(k,n3))

            for _ in [1,2]
                happyPhase1 = true

                if enableLog
                    println("ACA: tol = $(tol))")
                    println("grid $(n[1]) $(n[2]) $(n[3])")
                end
                
                # ACA 1
                t1 = time()
                T1e = evalTensor(Vector(1:n[1]),J,K,ff,vectorize) :: Array{Float64,3}
                t1 = time() - t1
                cf3F.numEvals += length(T1e)
                T1 = reshape(T1e,(n[1],r[2]*r[3])) :: Matrix{Float64}
                t2 = time()
                ~, tol = getTol(T1, pseudoLevel, tol)
                t2 = time() - t2
                t3 = time()
                Uc, ~, ~, I, I2 = ACA(T1, tol, n[1])
                t3 = time() - t3
                r[1] = size(I,1)
                JT1 = copy(J)
                KT1 = copy(K)

                # ACA 2
                T2e = evalTensor(I,Vector(1:n[2]),K,ff,vectorize) :: Array{Float64,3 }
                cf3F.numEvals += length(T2e) 
                T2 = reshape(permutedims(T2e,[2,1,3]),n[2],r[1]*r[3]) :: Matrix{Float64}
                ~, tol = getTol(T2, pseudoLevel, tol)
                Vc, ~, ~, J, J2 = ACA(T2, tol, n[2])
                r[2] = size(J,1)
                KT2 = copy(K)
                IT2 = copy(I)

                # ACA 3
                T3e = evalTensor(I,J,Vector(1:n[3]), ff,vectorize) :: Array{Float64,3}
                cf3F.numEvals += length(T3e)
                T3 = reshape(permutedims(T3e,[3,1,2]),n[3],r[1]*r[2]) :: Matrix{Float64}
                reltol, tol = getTol(T3, pseudoLevel, tol)
                Wc, ~, ~, K, K2 = ACA(T3, tol, n[3])
                r[3] = size(K,1)
                IT3 = copy(I)
                JT3 = copy(J)

                if enableLog
                    println("rank $(r[1]) $(r[2]) $(r[3]) evals $(cf3F.numEvals - evalsold)")
                    evalsold = cf3F.numEvals;
                end
                
                # Refine n
                breakFlag = false :: Bool
                for i in 1:3
                    while r[i] * 2 * sqrt(2) > n[i]
                        n[i] = reffun(n[i])
                        breakFlag = true
                    end
                end

                # Reinitialize after refinement
                if breakFlag
                    happyPhase1 = false;
                    if enableLog
                        println(" >> refine coarse grid")
                    end
                    break
                    # Proceed to refinement if r < 2
                elseif minimum(r) < 2
                    break
                end
            end 
            happyPhase1 = true
        end
        happy = true
    end
    return t1, t2, t3
end