using QuantumOptics
using LinearAlgebra
using Printf


const lib_ext = Sys.iswindows() ? ".dll" : Sys.isapple() ? ".dylib" : ".so"
const lib_path = joinpath(@__DIR__, "libmemory" * lib_ext)
function get_peak_memory()
    bytes = ccall((:get_peak_memory_bytes, lib_path), Clonglong, ())
    if bytes != -1
        @printf("Peak memory usage: %.4f GB\n", bytes / (1024^3))
    else
        println("Failed to retrieve memory usage.")
    end
end


function main(;mu_R_val::Float64 = 3.0, lambda_cav_val::Float64 = 0.1, set_name::String="SetII")
    start_time = time()
    get_peak_memory()

    delta = 1.0
    eps1 = delta
    eps2 = -delta
    en_gap = abs(eps1 - eps2)
    omega_t = 0.1 * en_gap
    omega_c = 0.2 * en_gap
    kappa1 = -0.2 * en_gap
    kappa2 = 0.3 * en_gap
    gamma_c = 0.05 * en_gap

    N_cav = 3
    frequency_factor = set_name == "SetII" ? 1.0/3 : 1.0
    omega_cav = en_gap * frequency_factor

    T = set_name == "SetI" ? 1.0/5 : 1.0
    mu_R = mu_R_val
    mu_L = -mu_R
    gammaL = 1.0
    gammaR = 1.0
    W = 10.0 * gammaL

    N1 = 20
    N2 = 10
    tau = 2.0

    println("Using $(N1) levels for vibrational modes")
    println("Using $(N_cav) levels for cavity modes")

    b_f1 = FockBasis(1)
    b_f2 = FockBasis(1)
    b_q1 = FockBasis(N1)
    b_q2 = FockBasis(N2)
    b_cav = FockBasis(N_cav)
    bases = [b_f1, b_f2, b_q1, b_q2, b_cav]
    b_tot = tensor(bases...)

    get_op(idx, op) = tensor([i == idx ? op : identityoperator(bases[i]) for i in 1:length(bases)]...)

    c1_sub = destroy(b_f1)
    c2_sub = destroy(b_f2)
    n1_sub = dagger(c1_sub) * c1_sub

    c1 = get_op(1, c1_sub)
    c2 = tensor(identityoperator(b_f1) - 2*n1_sub, c2_sub,
        identityoperator(b_q1), identityoperator(b_q2),
        identityoperator(b_cav))

    b1 = get_op(3, destroy(b_q1))
    b2 = get_op(4, destroy(b_q2))
    a_cav = get_op(5, destroy(b_cav))

    n1, n2 = dagger(c1)*c1, dagger(c2)*c2
    nb1, nb2 = dagger(b1)*b1, dagger(b2)*b2
    ncav = dagger(a_cav)*a_cav
    q1, q2 = (b1 + dagger(b1))/sqrt(2), (b2 + dagger(b2))/sqrt(2)

    #hamiltonian
    h0 = omega_t * (nb1 + 0.5*identityoperator(b_tot)) + omega_c * (nb2 + 0.5*identityoperator(b_tot))
    H_sys = (h0 + eps1*identityoperator(b_tot) + kappa1*q1) * n1 +
            (h0 + eps2*identityoperator(b_tot) + kappa2*q1) * n2 +
            gamma_c*q2*(dagger(c1)*c2 + dagger(c2)*c1) +
            omega_cav*ncav + lambda_cav_val*(dagger(c2)*c1 + dagger(c1)*c2)*(dagger(a_cav) + a_cav)

    f_FD(E, mu, temp) = 1.0 / (exp((E - mu) / temp) + 1.0)
    # J_spectral(E, mu) = (gammaL * W^2) / ((E - mu)^2 + W^2)
    fL1, fR2 = f_FD(eps1, mu_L, T), f_FD(eps2, mu_R, T)

    J = [dagger(c1), c1, dagger(c2), c2]
    Jdag = [dagger(j) for j in J]
    base_rates = [gammaL * fL1, gammaL * (1.0 - fL1), gammaR * fR2, gammaR * (1.0 - fR2)]

    f_dyn(t, rho) = (H_sys, J, Jdag, base_rates .* (1.0 - exp(-t/tau)))

    exp_ops = [n1, n2, nb1, nb2, ncav, c2*dagger(c2), c1*dagger(c1)]
    f_out(t, rho) = [real(expect(op, rho)) for op in exp_ops]

    rho0 = dm(tensor(fockstate(b_f1,0), fockstate(b_f2,0), fockstate(b_q1,0), fockstate(b_q2,0), fockstate(b_cav,0)))
    tspan = 0.0:0.2:20.0
    println("Running time propagation")
    tout, data = timeevolution.master_dynamic(tspan, rho0, f_dyn; fout=f_out)

    sw = 1.0 .- exp.(.-tout ./ tau)
    I_L = (sw .* (gammaL * fL1) .* [d[7] for d in data]) .- (sw .* (gammaL*(1-fL1)) .* [d[1] for d in data])
    I_R = (sw .* (gammaR * fR2) .* [d[6] for d in data]) .- (sw .* (gammaR*(1-fR2)) .* [d[2] for d in data])
    I_tot = (I_L .+ I_R) ./ 2.0


    pop1 = [d[1] for d in data]
    pop2 = [d[2] for d in data]
    nb1  = [d[3] for d in data]
    nb2  = [d[4] for d in data]
    ncav = [d[5] for d in data]
    current = I_tot


    current_file = joinpath(pwd(),
        @sprintf("lamda_%.4f", lambda_cav_val),
        @sprintf("mu_%.1f", mu_R), "current.masterEq.txt")
    mkpath(dirname(current_file))

    f = open(current_file, "w")
    for (i, val) in enumerate(tout)
        @printf(f, "%12.6f  %12.6f\n", val, I_tot[i])
    end
    close(f)

    population_file = joinpath(pwd(),
        @sprintf("lamda_%.4f", lambda_cav_val),
        @sprintf("mu_%.1f", mu_R), "population.masterEq.txt")
    mkpath(dirname(population_file))

    f = open(population_file, "w")
    for (i, val) in enumerate(tout)
        @printf(f, "%12.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f\n",
            val, pop1[i], pop2[i], nb1[i], nb2[i], ncav[i])
    end
    close(f)


    get_peak_memory()
    @printf("It took %.5f seconds\n", time()-start_time)
end

muR = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 3.0
lamda_cav = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.1

main(mu_R_val = muR, lambda_cav_val = lamda_cav)
