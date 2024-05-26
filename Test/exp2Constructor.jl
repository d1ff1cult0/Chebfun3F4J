using LinearAlgebra
using FastChebInterp
using InteractiveUtils
using Statistics

# Include helper functions
folder_path = joinpath(@__DIR__, "private")

files = readdir(folder_path)
for file in files
    include(joinpath(folder_path, file))
end

# Include data structures and related functions
include("chebfun3f.jl")

@time function cf3F2(f :: Function, enableLog::Bool = false, enableTime:: Bool = false, grid::Int64 = 17)
    # Informatie uit chebfunpref() van matlab
    pseudoLevel = eps()
    passSampleTest = true
    maxRestarts = 10

    if enableLog == 1
        println("\nChebfun3F constructor:")
        evalsold = 0;
    end

    #Initialize
    vectorize = false
    n = [grid,grid,grid]
    m = copy(n)
    r = [6,6,6]
    tol = pseudoLevel
    reltol = 0
    # Setup cf3F
    cf3F = chebfun3F([[]], [[]], [[]], [[]], [], 0, 0)
    cheb(i,n) = 0 .- (cos.((i.-1).*(pi/(n-1))))
    reffun(n) = (floor.(Int64,sqrt(2).^(floor.(2*log2.(n)) .+ 1)) .+ 1)

    # Check fast vectorization
    try
        A = f(1:2,1:2,1:2)
        vectorize = true
        if isScalar(A)
            g(x,y,z) = f(x,y,z) + 0*x + 0*y + 0*z # ensure vectorization
            f = g
        end
    catch
        vectorize = false
    end

    data = []

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
        happyPhase1 = false
        while !happyPhase1
            t1 = time()
            J = initializeIndexRandomly(r[2], n[2])
            K = initializeIndexRandomly(r[3], n[3])
            # Handle to evaluate tensor entries of T_c
            n1 = n[1]
            n2 = n[2]
            n3 = n[3]
            ff(i,j,k) = f(cheb(i,n1), cheb(j,n2), cheb(k,n3))

            for iterations in [1,2]
                happyPhase1 = true

                if enableLog
                    println("ACA: tol = $(tol))")
                    println("grid $(n[1]) $(n[2]) $(n[3])")
                end

                # ACA 1
                T1 = evalTensor(Vector(1:n[1]),J,K,ff,vectorize)
                cf3F.numEvals += numel(T1)
                T1 = reshape(T1,(n[1],r[2]*r[3]))
                ~, tol = getTol(T1, pseudoLevel, tol)
                Uc, ~, ~, I, I2 = ACA(T1, tol, n[1])
                r[1] = size(I,1)
                JT1 = copy(J)
                KT1 = copy(K)

                # ACA 2
                T2 = evalTensor(I,Vector(1:n[2]),K,ff,vectorize)
                cf3F.numEvals += numel(T2) 
                T2 = reshape(permutedims(T2,[2,1,3]),n[2],r[1]*r[3])
                ~, tol = getTol(T2, pseudoLevel, tol)
                Vc, ~, ~, J, J2 = ACA(T2, tol, n[2])
                r[2] = size(J,1)
                KT2 = copy(K)
                IT2 = copy(I)

                # ACA 3
                T3 = evalTensor(I,J,Vector(1:n[3]), ff,vectorize)
                cf3F.numEvals += numel(T3)
                T3 = reshape(permutedims(T3,[3,1,2]),n[3],r[1]*r[2])
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
                breakFlag = false
                while r[1]*2*sqrt(2) > n[1]
                    n[1] = reffun(n[1])
                    breakFlag = true
                end
                while r[2]*2*sqrt(2) > n[2]
                    n[2] = reffun(n[2])
                    for i = 1:3
                        while r[i]*2*sqrt(2) > n[i]
                            n[i] = reffun(n[i])
                        end
                    end
                    breakFlag = true
                end
                while r[3]*2*sqrt(2) > n[3]
                    n[3] = reffun(n[3])
                    breakFlag = true
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
            if enableTime
                println("Phase 1 time: $(time() - t1)")
            end
        end

        # Phase 2
        # Catch the rank zero function:
        t2 = time()
        if size(I,1) == 0 || size(J,1) == 0 || size(K,1) == 0
            error("rank zero matrix")
        else
            # Refine
            m = copy(n)

            Uf = Uc
            ct2 = createCT2(Uf)
            resolvedU = happinessCheck(ct2, tol)
            if !resolvedU
                m[1] = 2*m[1]-1
            end

            Vf = Vc
            ct2 = createCT2(Vf)
            resolvedV = happinessCheck(ct2, tol)
            if !resolvedV
                m[2] = 2*m[2]-1
            end

            Wf = Wc
            ct2 = createCT2(Wf)
            resolvedW = happinessCheck(ct2, tol)
            if !resolvedW
                m[3] = 2*m[3]-1
            end

            # Add function evaluations and check again
            refiter = 0
             while !resolvedU || !resolvedV || !resolvedW
                refiter += 1
                # function handle to evaluate T_f
                m1 = m[1]
                m2 = m[2]
                m3 = m[3]
                ff(i,j,k) = f(cheb(i,m1),cheb(j,m2),cheb(k,m3))

                # Map the indices from T_c to T_f       
                refFactor = [0, 0, 0]
                iter = 0
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
                Jr = ref(JT1, refFactor[2])
                Kr = ref(KT1, refFactor[3])
                if !resolvedU
                    Uold = deepcopy(Uf)
                    Ij, Ik = ind2sub(size(Jr, 1), I2)
                    if vectorize
                        n1 = Int(floor(m[1]/2)) ; n2 = size(I2,1)
                        x = hcat(2:2:m[1])
                        X = repeat(x, 1, n2)
                        y = Jr[Ij]'
                        Y = repeat(y, n1, 1)
                        z = Kr[Ik]'
                        Z = repeat(z, n1, 1)
                        Uf = vcat(Uf, zeros(Float64, m[1] - size(Uf, 1), size(Uf,2)))
                        Uf[2:2:m[1], 1:n2] .= ff(X, Y, Z)
                        Uf[1:2:m[1], 1:n2] .= Uold
                        cf3F.numEvals += numel(X)
                    else
                        Uf = zeros(m[1], size(I2, 2))
                        for i in 1:m[1]
                            for j in 1:size(I2, 2)
                                if i % 2 == 0
                                    Uf[i, j] = ff(i, Jr[Ij[j]], Kr[Ik[j]])
                                    cf3F.numEvals += 1
                                else
                                    Uf[i, j] = Uold[(i + 1) / 2, j]
                                end
                            end
                        end
                    end
                    ct2 = createCT2(Uf)
                    resolvedU = happinessCheck(ct2, reltol)
                    if !resolvedU
                        m[1] = 2 * m[1] - 1
                    end
                end

                # V
                Ir = ref(IT2, refFactor[1])
                Kr = ref(KT2, refFactor[3])
                if !resolvedV
                    Vold = deepcopy(Vf)
                    Ji, Jk = ind2sub(size(Ir, 1), J2)
                    if vectorize
                        n1 = Int(floor(m[2]/2)) ; n2 = size(J2,1)
                        x = Ir[Ji]'
                        X = repeat(x, n1, 1)
                        y = hcat(2:2:m[2])
                        Y = repeat(y, 1, n2)
                        z = Kr[Jk]'
                        Z = repeat(z, n1, 1)
                        Vf = vcat(Vf, zeros(Float64, m[2] - size(Vf, 1), size(Vf,2)))
                        Vf[2:2:m[2], 1:n2] .= ff(X, Y, Z)
                        Vf[1:2:m[2], 1:n2] .= Vold
                        cf3F.numEvals += numel(X)
                    else
                        Vf = zeros(m[2], size(J2, 2))
                        for i in 1:m[2]
                            for j in 1:size(J2, 2)
                                if i % 2 == 0
                                    Vf[i, j] = ff(Ir[Ji[j]], i, Kr[Jk[j]])
                                    cf3F.numEvals += 1
                                else
                                    Vf[i, j] = Vold[(i+1) / 2, j]
                                end
                            end
                        end
                    end
                    ct2 = createCT2(Vf)
                    resolvedV = happinessCheck(ct2, reltol)
                    if !resolvedV
                        m[2] = 2 * m[2] - 1
                    end
                end

                # W
                Ir = ref(IT3, refFactor[1])
                Jr = ref(JT3, refFactor[2])
                if !resolvedW
                    Wold = deepcopy(Wf)
                    Ki, Kj = ind2sub(size(Ir, 1), K2)
                    if vectorize
                        n1 = Int(floor(m[3]/2)) ; n2 = size(K2,1)
                        x = Ir[Ki]'
                        X = repeat(x, n1, 1)
                        y = Jr[Kj]'
                        Y = repeat(y, n1, 1)
                        z = hcat(2:2:m[3])
                        Z = repeat(z, 1, n2)
                        Wf = vcat(Wf, zeros(Float64, m[3] - size(Wf, 1), size(Wf,2)))
                        Wf[2:2:m[3], 1:n2] .= ff(X, Y, Z)
                        Wf[1:2:m[3], 1:n2] .= Wold
                        cf3F.numEvals += numel(X)
                    else
                        Wf = zeros(m[3], size(K2, 2))
                        for i in 1:m[3]
                            for j in 1:size(K2, 2)
                                if i % 2 == 0
                                    Wf[i, j] = ff(Ir[Ki[j]], Jr[Kj[j]], i)
                                    cf3F.numEvals += 1
                                else
                                    Wf[i, j] = Wold[(i+1)÷2, j]
                                end
                            end
                        end
                    end
                    ct2 = createCT2(Wf)
                    resolvedW = happinessCheck(ct2, reltol)
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

            if enableTime
                println("Phase 2 time: $(time() - t2)")
            end

            # Phase 3
            t3 = time()
            Q1 = Matrix(qr(Uf).Q)
            I = DEIM(Q1)
            U = Q1/Q1[I,:]

            Q2 = Matrix(qr(Vf).Q)
            J = DEIM(Q2)
            V = Q2/Q2[J,:]

            Q3 = Matrix(qr(Wf).Q)
            K = DEIM(Q3)
            W = Q3/Q3[K,:]

            cf3F.r = [size(U, 2), size(V, 2), size(W, 2)]

            cf3F.U = chebfun(U)
            cf3F.V = chebfun(V)
            cf3F.W = chebfun(W)

            m1 = m[1]
            m2 = m[2]
            m3 = m[3]
            ff2(i,j,k) = f(cheb(i,m1),cheb(j,m2),cheb(k,m3))

            cf3F.C = evalTensor(I,J,K,ff2,vectorize)
            cf3F.numEvals += numel(cf3F.C)

            if enableLog == 1
                println("Core & Assembling: max cond $(maximum([cond(Q1[I, :]), cond(Q2[J, :]), cond(Q3[K, :])])) evals $(cf3F.numEvals-evalsold)")
            end
            evalsold = cf3F.numEvals
            if enableTime
                println("Phase 3 time: $(time() - t3)")
            end
        end

        # Sample Test
        if passSampleTest && cf3F.numRestarts < maxRestarts
            happy, maxError, meanError = sampleTest(cf3F, f, tol, enableLog)
            push!(data,(n,maxError,meanError))
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

    return cf3F, data
end

