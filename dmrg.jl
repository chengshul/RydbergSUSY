using ITensors, HDF5

function dmrg0(L::Int, r1::Float64, r2::Float64, Δ::Float64, C0::Float64,
        C1::Float64, C2::Float64, chi::Int)
    f1 = h5open("../../../data/data0/psi0_L_$(L)_r1_$(r1).h5","r")
    psi0_init = read(f1,"psi",MPS)
    close(f1)
    chi0 = 60

    sites = siteinds(psi0_init)
    ampo = OpSum()
    #Rb = 1.0
    L0 = L÷2
    for j in 0:L0-1
        j1 = 4*j+1
        j2 = 4*j+2
        j3 = 4*j+3
        j4 = 4*j+4
        j5 = 4*((j+1)%L0)+1
        j6 = 4*((j+1)%L0)+2
        j7 = 4*((j+1)%L0)+3
        j8 = 4*((j+1)%L0)+4
        j9 = 4*((j+2)%L0)+1
        j10 = 4*((j+2)%L0)+2

        ampo += 1/2, "X",j1
        ampo += 1/2, "X",j2
        ampo += 1/2, "X",j3
        ampo += 1/2, "X",j4

        ampo += C1/r2^6, "ProjUp",j1,"ProjUp",j2
        ampo += C2/r2^6, "ProjUp",j3,"ProjUp",j4
        ampo += C0/r1^6, "ProjUp",j1,"ProjUp",j3
        ampo += C0/r1^6, "ProjUp",j2,"ProjUp",j4
        ampo += C0/r1^6, "ProjUp",j3,"ProjUp",j5
        ampo += C0/r1^6, "ProjUp",j4,"ProjUp",j6

        ampo += C0/(r1^2+r2^2)^3, "ProjUp",j1,"ProjUp",j4
        ampo += C0/(r1^2+r2^2)^3, "ProjUp",j2,"ProjUp",j3
        ampo += C0/(r1^2+r2^2)^3, "ProjUp",j3,"ProjUp",j6
        ampo += C0/(r1^2+r2^2)^3, "ProjUp",j4,"ProjUp",j5
        
        ampo += C1/(2*r1)^6, "ProjUp",j1,"ProjUp",j5
        ampo += C1/(2*r1)^6, "ProjUp",j2,"ProjUp",j6
        ampo += C2/(2*r1)^6, "ProjUp",j3,"ProjUp",j7
        ampo += C2/(2*r1)^6, "ProjUp",j4,"ProjUp",j8
        
        ampo += C1/(4*r1^2+r2^2)^3, "ProjUp",j1,"ProjUp",j6
        ampo += C1/(4*r1^2+r2^2)^3, "ProjUp",j2,"ProjUp",j5
        ampo += C2/(4*r1^2+r2^2)^3, "ProjUp",j3,"ProjUp",j8
        ampo += C2/(4*r1^2+r2^2)^3, "ProjUp",j4,"ProjUp",j7

        ampo += C0/(3*r1)^6, "ProjUp",j1,"ProjUp",j7
        ampo += C0/(3*r1)^6, "ProjUp",j2,"ProjUp",j8
        ampo += C0/(3*r1)^6, "ProjUp",j3,"ProjUp",j9
        ampo += C0/(3*r1)^6, "ProjUp",j4,"ProjUp",j10

        ampo += C0/(9*r1^2+r2^2)^3, "ProjUp",j1,"ProjUp",j8
        ampo += C0/(9*r1^2+r2^2)^3, "ProjUp",j2,"ProjUp",j7
        ampo += C0/(9*r1^2+r2^2)^3, "ProjUp",j3,"ProjUp",j10
        ampo += C0/(9*r1^2+r2^2)^3, "ProjUp",j4,"ProjUp",j9

        ampo += Δ, "ProjUp",j1
        ampo += Δ, "ProjUp",j2
        ampo += Δ, "ProjUp",j3
        ampo += Δ, "ProjUp",j4
    end
    ampo += -0.0001,"ProjDn",1
    H = MPO(ampo,sites)

    sweeps = Sweeps(100)
    setnoise!(sweeps, 1e-6,1e-6,1e-7,1e-7,1e-8,1e-8,0)
    setmaxdim!(sweeps, chi)
    setmindim!(sweeps, chi0)
    setcutoff!(sweeps, 1E-12)

    E0, psi0 = dmrg(H,psi0_init, sweeps)
    println("Final energy = $(E0)")

    f1 = h5open("../../../data/data0/H_2a_L_$(L)_Delta_$(Δ)_r1_$(r1)_r2_$(r2).h5","w")
    write(f1,"H",H)
    close(f1)

    f2 = h5open("../../../data/data0/psi_2a_L_$(L)_Delta_$(Δ)_r1_$(r1)_r2_$(r2)_chi_$(chi).h5","w")
    write(f2,"psi",psi0)
    close(f2)

    f2 = h5open("../../../data/data0/psi0_L_$(L)_r1_$(r1).h5","w")
    write(f2,"psi",psi0)
    close(f2)

    magz = expect(psi0,"Z")
    corr = correlation_matrix(psi0,"Z","Z")
    corr_f = zeros(L-1,4)
    for i in 1:L-1
        op_i0 = [[1,2*i+1],[2,2*i+1],[1,2*i+2],[2,2*i+2]]
        for i0 in 1:4
            os = Tuple{String, Int64, Vararg{Int64}}[]
            push!(os,("Z",op_i0[i0][1]))
            for j in 1:i-1
                push!(os,("Swap",2*j+1,2*j+2))
            end
            push!(os,("Z",op_i0[i0][2]))
            gates = ops(os,sites)
            ψ = apply(gates, psi0; cutoff = 1e-15)
            corr_f[i,i0] = inner(psi0,ψ)
        end
    end    
    
    M_0 = OpSum()
    for j in 0:L-1
        j1 = 2*j+1
        j2 = 2*j+2
        M_0 += 1, "Z",j1
        M_0 += -1, "Z",j2
    end
    M = MPO(M_0,sites)
    psi1 = apply(M,psi0)
    psi1 = apply(M,psi1)
    psi1 = apply(M,psi1)
    psi1 = apply(M,psi1)
    M4 = inner(psi0,psi1)

    SvN_l = zeros(2*L-1)

    orthogonalize!(psi0, 1)
    U,S,V = svd(psi0[1], (siteind(psi0,1)))
    SvN = 0.0
    for n in 1:dim(S, 1)
        p = S[n,n]^2
        SvN -= p * log(p)
    end
    SvN_l[1] = SvN

    for b in 2:2*L-1 # = N÷2
        orthogonalize!(psi0, b)
        U,S,V = svd(psi0[b], (linkind(psi0, b-1), siteind(psi0,b)))
        SvN = 0.0
        for n in 1:dim(S, 1)
            p = S[n,n]^2
            SvN -= p * log(p)
        end
        SvN_l[b] = SvN
    end

    return magz, corr, corr_f, SvN_l, M4
end

function main(L::Int, r1::Float64, r2::Float64, Δ1::Float64, Δ2::Float64, chi::Int)
    N = 2*L
    N_Δ = 21
    C1 = 1.0
    C2 = 1.64
    C0 = -0.97
    #r2 = 0.8
    #N_Rb = 2

    #E_all = zeros(Float64,N_Δ,N_Rb)
    corr1_all = zeros(Float64,N_Δ,N)
    corr2_all = zeros(Float64,N_Δ,N,N)
    corrf_all = zeros(Float64,N_Δ,L-1,4)
    SvN_all = zeros(Float64,N_Δ,N-1)
    M4_all = zeros(Float64,N_Δ)

    Δ_l = range(Δ1,Δ2; length=N_Δ)
    #Rb_l = range(2,4; length=N_Rb)

    for i_Δ in 1:N_Δ
        Δ = Δ_l[i_Δ]
        magz, corr, corrf, SvN, M4 = dmrg0(L, r1, r2, Δ, C0, C1, C2, chi)
        corr1_all[i_Δ,:] = magz
        corr2_all[i_Δ,:,:] = corr
        corrf_all[i_Δ,:,:] = corrf
        SvN_all[i_Δ,:,:] = SvN
        M4_all[i_Δ] = M4
    end
    return corr1_all, corr2_all, corrf_all, SvN_all, M4_all
end

L = parse(Int, ARGS[1])
r1 = parse(Float64, ARGS[2])
r2 = 0.5
Δ1 = parse(Float64, ARGS[3])
Δ2 = parse(Float64, ARGS[4])
chi = parse(Int, ARGS[5])

corr1_all, corr2_all, corrf_all, SvN_all, M4_all = main(L, r1, r2, Δ1, Δ2, chi)

open("../../../data/corr1_2a_L_$(L)_r1_$(r1)_r2_$(r2)_D1_$(Δ1)_D2_$(Δ2)_chi_$(chi).dat", "w") do io
    write(io, corr1_all)
end
open("../../../data/corr2_2a_L_$(L)_r1_$(r1)_r2_$(r2)_D1_$(Δ1)_D2_$(Δ2)_chi_$(chi).dat", "w") do io
    write(io, corr2_all)
end
open("../../../data/corrf_2a_L_$(L)_r1_$(r1)_r2_$(r2)_D1_$(Δ1)_D2_$(Δ2)_chi_$(chi).dat", "w") do io
    write(io, corrf_all)
end
open("../../../data/S_2a_L_$(L)_r1_$(r1)_r2_$(r2)_D1_$(Δ1)_D2_$(Δ2)_chi_$(chi).dat", "w") do io
    write(io, SvN_all)
end
open("../../../data/M4_2a_L_$(L)_r1_$(r1)_r2_$(r2)_D1_$(Δ1)_D2_$(Δ2)_chi_$(chi).dat", "w") do io
    write(io, M4_all)
end

