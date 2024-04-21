using ITensors, HDF5

function dmrg0(L::Int, r1::Float64, r2::Float64, Δ::Float64, C0::Float64,
        C1::Float64, C2::Float64)
    sites = siteinds("S=1/2",2*L) #; conserve_qns=true)
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

    #psi0_init = randomMPS(sites)
    state = ["↑" for n=1:2*L]
    for i in 1:2:2*L
        state[i] = "↓"
    end
    psi0_init = productMPS(sites,state)

    sweeps = Sweeps(200)
    setnoise!(sweeps, 1e-6,1e-6,1e-7,1e-7,1e-8,1e-8,0)
    setmaxdim!(sweeps, 10,10,20,20,50,50,50,100,100,100,
                       200)
    setcutoff!(sweeps, 1E-12)
    chi = sweeps.maxdim[end]

    E0, psi0 = dmrg(H,psi0_init, sweeps)
    println("Final energy = $(E0)")

    f2 = h5open("../../../data/data0/psi0_L_$(L)_r1_$(r1).h5","w")
    write(f2,"psi",psi0)
    close(f2)

    return
end

L = parse(Int, ARGS[1])
r1 = parse(Float64, ARGS[2])
r2 = 0.5
Δ = parse(Float64, ARGS[3])

C1 = 1.0
C2 = 1.64
C0 = -0.97
dmrg0(L,r1,r2,Δ,C0,C1,C2)


