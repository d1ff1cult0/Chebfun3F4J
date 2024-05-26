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

@time function cf3F(f :: Function, enableLog::Bool = false, enableTime:: Bool = false, grid::Int64 = 17)
    # Informatie uit chebfunpref() van matlab
    pseudoLevel = eps() :: Float64
    passSampleTest = true :: Bool
    maxRestarts = 10 :: Int64

    if enableLog == 1
        println("\nChebfun3F constructor:")
        evalsold = 0;
    end

    #Initialize
    vectorize = false :: Bool
    n = fill(grid,3) :: Vector{Int64}
    m = copy(n) :: Vector{Int64}
    r = fill(6,3) :: Vector{Int64}
    tol = pseudoLevel :: Float64
    reltol = 0.0 :: Float64
    # Setup cf3F
    cf3F = chebfun3F([[]], [[]], [[]], [[]], [], 0, 0) :: chebfun3F
    cheb(i,n) = 0 .- (cos.((i.-1).*(pi/(n-1))))
    reffun(n) = (floor.(Int64,sqrt(2).^(floor.(2*log2.(n)) .+ 1)) .+ 1)

    t1_times = []
    t2_times = []
    t3_times = []

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

    #Main loop
    happy::Bool = false
    while !happy
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
            t1 = time()
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
                T1e = evalTensor(Vector(1:n[1]),J,K,ff,vectorize) :: Array{Float64,3}
                cf3F.numEvals += length(T1e)
                T1 = reshape(T1e,(n[1],r[2]*r[3])) :: Matrix{Float64}
                ~, tol = getTol(T1, pseudoLevel, tol)
                Uc, ~, ~, I, I2 = ACA(T1, tol, n[1])
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
            t1_times = push!(t1_times, time() - t1)
            if enableTime
                println("Phase 1 time: $(t1_times[end])")
            end
        end

        # Phase 2
        # Catch the rank zero function:
        t2 = time()
        if size(I,1) == 0 || size(J,1) == 0 || size(K,1) == 0
            error("rank zero matrix")
        else
            # Refine
            m = copy(n) :: Vector{Int64}

            Uf = Uc :: Matrix{Float64}
            Vf = Vc :: Matrix{Float64}
            Wf = Wc :: Matrix{Float64}

            resolvedU = happinessCheck(createCT2(Uf), tol) :: Bool
            if !resolvedU
                m[1] = 2*m[1]-1
            end
            
            resolvedV = happinessCheck(createCT2(Vf), tol) :: Bool
            if !resolvedV
                m[2] = 2*m[2]-1
            end
            
            resolvedW = happinessCheck(createCT2(Wf), tol) :: Bool
            if !resolvedW
                m[3] = 2*m[3]-1
            end

            # Add function evaluations and check again
            refiter = 0 :: Int64
            while !(resolvedU && resolvedV && resolvedW)
                refiter += 1
                # function handle to evaluate T_f
                m1 :: Int64, m2 :: Int64, m3 :: Int64 = m
                ff(i,j,k) = f(cheb(i,m1),cheb(j,m2),cheb(k,m3))

                # Map the indices from T_c to T_f       
                refFactor = zeros(Int64,3) :: Vector{Int64}
                iter = 0 :: Int64
                while minimum(refFactor) == 0
                    iter += 1
                    if n[1] * iter - (iter - 1) == m[1]
                        refFactor[1] = iter
                    end
                    if n[2] * iter - (iter - 1) == m[2]
                        refFactor[2] = iter
                    end
                    if n[3] * iter - (iter - 1) == m[3]
                        refFactor[3] = iter
                    end
                end
                ref(i, r) = r .* i .- (r .- 1)
                
                # U
                Jr = ref(JT1, refFactor[2]) :: Vector{Int64}
                Kr = ref(KT1, refFactor[3]) :: Vector{Int64}
                if !resolvedU
                    Uold = deepcopy(Uf) :: Matrix{Float64}
                    Ij :: Vector{Int64}, Ik :: Vector{Int64} = ind2sub(size(Jr, 1), I2)
                    if vectorize
                        n1 = div(m[1], 2)
                        n2 = size(I2, 1)
                        x = collect(2:2:m[1]) :: Vector{Int64}
                        y = Jr[Ij] :: Vector{Int64}
                        z = Kr[Ik] :: Vector{Int64}
                        
                        Uf = vcat(Uf, zeros(Float64, m[1] - size(Uf, 1), size(Uf, 2)))
                        Uf[2:2:m[1], :] .= ff.(x, repeat(y', n1, 1), repeat(z', n1, 1))
                        Uf[1:2:m[1], :] .= Uold
                        cf3F.numEvals += length(x) * length(y) * length(z)
                    else
                        newUf = zeros(Float64, m[1], size(Uold, 2))
                        for i in 1:m[1]
                            for j in 1:size(I2, 2)
                                if iseven(i)
                                    newUf[i, j] = ff(i, Jr[Ij[j]], Kr[Ik[j]])
                                    cf3F.numEvals += 1
                                else
                                    newUf[i, j] = Uold[div(i + 1, 2), j]
                                end
                            end
                        end
                        Uf = newUf
                    end
                    resolvedU = happinessCheck(createCT2(Uf), reltol)
                    if !resolvedU
                        m[1] = 2 * m[1] - 1
                    end
                end

                # V
                Ir = ref(IT2, refFactor[1]) :: Vector{Int64}
                Kr = ref(KT2, refFactor[3]) :: Vector{Int64}
                if !resolvedV
                    Vold = deepcopy(Vf) :: Matrix{Float64}
                    Ji, Jk = ind2sub(size(Ir, 1), J2) :: Tuple{Vector{Int64}, Vector{Int64}}
                    
                    if vectorize
                        n1 = div(m[2], 2)
                        n2 = size(J2, 1)
                        x = Ir[Ji] :: Vector{Int64}
                        y = collect(2:2:m[2]) :: Vector{Int64}
                        z = Kr[Jk] :: Vector{Int64}
                        
                        Vf = vcat(Vf, zeros(Float64, m[2] - size(Vf, 1), size(Vf, 2)))
                        Vf[2:2:m[2], :] .= ff.(repeat(x', n1, 1), y, repeat(z', n1, 1))
                        Vf[1:2:m[2], :] .= Vold
                        cf3F.numEvals += length(x) * length(y) * length(z)
                    else
                        newVf = zeros(Float64, m[2], size(Vold, 2))
                        for i in 1:m[2]
                            for j in 1:size(J2, 2)
                                if iseven(i)
                                    Vf[i, j] = ff(Ir[Ji[j]], i, Kr[Jk[j]])
                                    cf3F.numEvals += 1
                                else
                                    Vf[i, j] = Vold[div(i + 1, 2), j]
                                end
                            end
                        end
                        Vf = newVf
                    end

                    resolvedV = happinessCheck(createCT2(Vf), reltol)
                    if !resolvedV
                        m[2] = 2 * m[2] - 1
                    end
                end

                # W
                Ir = ref(IT3, refFactor[1]) :: Vector{Int64}
                Jr = ref(JT3, refFactor[2]) :: Vector{Int64}
                if !resolvedW
                    Wold = deepcopy(Wf) :: Matrix{Float64}
                    Ki :: Vector{Int64}, Kj :: Vector{Int64} = ind2sub(size(Ir, 1), K2)
                    if vectorize
                        n1 = div(m[3], 2)
                        n2 = size(K2, 1)
                        x = Ir[Ki] :: Vector{Int64}
                        y = Jr[Kj] :: Vector{Int64}
                        z = collect(2:2:m[3]) :: Vector{Int64}
                        
                        Wf = vcat(Wf, zeros(Float64, m[3] - size(Wf, 1), size(Wf, 2)))
                        Wf[2:2:m[3], :] .= ff.(repeat(x', n1, 1), repeat(y', n1, 1), repeat(z, 1, n2))
                        Wf[1:2:m[3], :] .= Wold
                        cf3F.numEvals += length(x) * length(y) * length(z)
                    else
                        newWf = zeros(Float64, m[3], size(Wold, 2))
                        for i in 1:m[3]
                            for j in 1:size(K2, 2)
                                if iseven(i)
                                    Wf[i, j] = ff(Ir[Ki[j]], Jr[Kj[j]], i)
                                    cf3F.numEvals += 1
                                else
                                    Wf[i, j] = Wold[div(i + 1, 2), j]
                                end
                            end
                        end
                        Wf = newWf
                    end
                    resolvedW = happinessCheck(createCT2(Wf), reltol)
                    if !resolvedW
                        m[3] = 2 * m[3] - 1
                    end
                end
            end

            if enableLog
                println("Refinement: $(m[1]) $(m[2]) $(m[3]) evals $(cf3F.numEvals - evalsold)")
                evalsold = cf3F.numEvals;
            end
            ~, tol = getTol(Uf, pseudoLevel, tol);
            ~, tol = getTol(Vf, pseudoLevel, tol);
            ~, tol = getTol(Wf, pseudoLevel, tol);

            t2_times = push!(t2_times, time() - t2)
            if enableTime
                println("Phase 2 time: $(t2_times[end])")
            end

            # Phase 3
            t3 = time()
            Q1 = Matrix(qr(Uf).Q) :: Matrix{Float64}
            I = DEIM(Q1) :: Vector{Int64}
            U = Q1/Q1[I,:] :: Matrix{Float64}

            Q2 = Matrix(qr(Vf).Q) :: Matrix{Float64}
            J = DEIM(Q2) :: Vector{Int64}
            V = Q2/Q2[J,:] :: Matrix{Float64}

            Q3 = Matrix(qr(Wf).Q) :: Matrix{Float64}
            K = DEIM(Q3) :: Vector{Int64}
            W = Q3/Q3[K,:] :: Matrix{Float64}

            cf3F.r = [size(U, 2), size(V, 2), size(W, 2)]

            cf3F.U = chebfun(U)
            cf3F.V = chebfun(V)
            cf3F.W = chebfun(W)

            m1, m2, m3 = m
            ff2(i,j,k) = f(cheb(i,m1),cheb(j,m2),cheb(k,m3))

            cf3F.C = evalTensor(I,J,K,ff2,vectorize)
            cf3F.numEvals += length(cf3F.C)

            if enableLog == 1
                println("Core & Assembling: max cond $(maximum([cond(Q1[I, :]), cond(Q2[J, :]), cond(Q3[K, :])])) evals $(cf3F.numEvals-evalsold)")
            end
            evalsold = cf3F.numEvals
            
            t3_times = push!(t3_times, time() - t3)
            if enableTime
                println("Phase 3 time: $(t3_times[end])")
            end
        end

        # Sample Test
        if passSampleTest && cf3F.numRestarts < maxRestarts
            happy, ~, ~ = sampleTest(cf3F, f, tol, enableLog)
            cf3F.numEvals += 30
            if enableLog
                println("evals $(cf3F.numEvals - evalsold)")
            end
            evalsold = cf3F.numEvals
            if happy && enableLog
                println("total $(cf3F.numEvals)")
            end
        else
            happy = true
            if enableLog
                println("total $(cf3F.numEvals)")
            end
        end

        # Restart
        if !happy
            if cf3F.numRestarts + 1 == maxRestarts
                println("Warning: Chebfun3F: max number of restarts reached")
                return
            end

            # Increase n
            n = reffun(n)
            cf3F.numRestarts += 1

            if enableLog
                println(">> restart")
            end

            # Ensure r is large enough for
            # (1,r,r) functions
            if r[1] > 1 || r[2] > 1 || r[3] > 1
                if r[1] < 3
                    r[3] = max(6,2*r[3])
                    r[2] = max(6,2*r[2])
                elseif r[2] < 2
                    r[1] = max(6,2*r[1])
                    r[3] = max(6,2*r[3])
                elseif r[3] < 2
                    r[1] = max(6,2*r[1])
                    r[2] = max(6,2*r[2])
                end            
            end

            # Very low-rank functions
            r[1] = max(r[1],3)
            r[2] = max(r[2],3)
            r[3] = max(r[3],3)

            if cf3F.numRestarts >= maxRestarts/2
                r *= 2
            end
        end
    end
    if enableLog 
        println("succes")
    end

    println(t1_times)
    println(t2_times)
    println(t3_times)
    return cf3F, mean(t1_times), mean(t2_times), mean(t3_times)
end 

function isScalar(x)
    return isa(x, Number) || isa(x, String)
end

