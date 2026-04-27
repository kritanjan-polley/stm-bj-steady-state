using QuantumToolbox
using LinearAlgebra
using Printf

const lib_ext = Sys.iswindows() ? ".dll" : Sys.isapple() ? ".dylib" : ".so"
const lib_path = joinpath(@__DIR__, "libmemory" * lib_ext)
function get_peak_memory()
    try
        bytes = ccall((:get_peak_memory_bytes, lib_path), Clonglong, ())
        if bytes != -1
            @printf("Peak memory usage: %.4f GB\n", bytes / (1024^3))
        else
            println("Failed to retrieve memory usage.")
        end
    catch e
        println("Memory tracking library not found or failed. Skipping memory check.")
    end
end

function main(;mu_R_val::Float64 = 3.0, lambda_cav_val::Float64 = 0.1, set_name="SetII")
    start_time = time()
    get_peak_memory()

    delta = 1.0
    eps1 = delta
    eps2 = -delta
    en_gap = abs(eps1 - eps2)

    N_cav = 20
    frequency_factor = set_name == "SetII" ? 1.0/3 : 1.0
    omega_cav = en_gap * frequency_factor

    T = set_name == "SetI" ? 1.0/5 : set_name == "SetIV" ? 1.0/2 : 1.0
    mu_R = -mu_R_val
    mu_L = mu_R_val
    gammaL = 1.0
    gammaR = 1.0
    W = 10.0 * gammaL
    tau = 2.0

    println("Using $(N_cav) levels for cavity mode")
    println("Vibrational modes (tuning and coupling) have been removed.")

    println("="^50)
    @printf("Energy levels: %.4f and %.4f\n", eps1, eps2)
    @printf("Cavity mode frequency = %.4f\n", omega_cav)
    @printf("Cavity coupling lambda = %.4f\n", lambda_cav_val)
    @printf("Bias voltages are: %.4f and %.4f\n", mu_L, mu_R)
    @printf("temperature = %.4f\n", T)
    println("="^50)

    dims = [2, 2, N_cav]
    get_op(idx, op) = tensor([i == idx ? op : qeye(dims[i]) for i in 1:length(dims)]...)

    c1_sub = destroy(2)
    c2_sub = destroy(2)
    n1_sub = c1_sub' * c1_sub

    c1 = get_op(1, c1_sub)
    c2 = tensor(qeye(2) - 2*n1_sub, c2_sub, qeye(N_cav))

    a_cav = get_op(3, destroy(N_cav))

    n1, n2 = c1'*c1, c2'*c2
    ncav = a_cav'*a_cav

    H_sys = eps1 * n1 +
            eps2 * n2 +
            omega_cav * ncav +
            lambda_cav_val * (c2'*c1 + c1'*c2) * (a_cav' + a_cav)

    f_FD(E, mu, temp) = 1.0 / (exp((E - mu) / temp) + 1.0)
    fL1, fR2 = f_FD(eps1, mu_L, T), f_FD(eps2, mu_R, T)

    J = [c1', c1, c2', c2]
    base_rates = [gammaL * fL1, gammaL * (1.0 - fL1), gammaR * fR2, gammaR * (1.0 - fR2)]
    c_ops = [QobjEvo( (J[i], (p, t) -> sqrt(base_rates[i] * (1.0 - exp(-t/tau)))) ) for i in 1:4]

    exp_ops = [n1, n2, ncav, c2*c2', c1*c1']
    psi0 = tensor(fock(2,0), fock(2,0), fock(N_cav,0))
    rho0 = ket2dm(psi0)

    tspan = 0.0:0.2:20.0
    println("Running time propagation...")

    sol = mesolve(H_sys, rho0, tspan, c_ops, e_ops=exp_ops)

    d_n1    = real.(sol.expect[1, :])
    d_n2    = real.(sol.expect[2, :])
    d_ncav  = real.(sol.expect[3, :])
    d_hole2 = real.(sol.expect[4, :])
    d_hole1 = real.(sol.expect[5, :])

    sw = 1.0 .- exp.(.-tspan ./ tau)
    I_L = (sw .* (gammaL * fL1) .* d_hole1) .- (sw .* (gammaL*(1-fL1)) .* d_n1)
    I_R = (sw .* (gammaR * fR2) .* d_hole2) .- (sw .* (gammaR*(1-fR2)) .* d_n2)
    I_tot = (I_L .- I_R) ./ 2.0

    current_file = joinpath(pwd(), set_name,
        @sprintf("lamda_%.4f", lambda_cav_val),
        @sprintf("mu_%.1f", mu_L), "current.masterEq.txt")
    # current_file = joinpath(pwd(), "current.masterEq.txt")
    mkpath(dirname(current_file))

    f = open(current_file, "w")
    for (i, val) in enumerate(tspan)
        @printf(f, "%12.6f  %12.6f %12.6f\n", val, I_L[i], I_R[i])
    end
    close(f)

    population_file = joinpath(pwd(), set_name,
        @sprintf("lamda_%.4f", lambda_cav_val),
        @sprintf("mu_%.1f", mu_L), "population.masterEq.txt")
    # population_file = joinpath(pwd(), "population.masterEq.txt")
    mkpath(dirname(population_file))

    f = open(population_file, "w")
    for (i, val) in enumerate(tspan)
        @printf(f, "%12.6f  %12.6f  %12.6f  %12.6f\n",
            val, d_n1[i], d_n2[i], d_ncav[i])
    end
    close(f)

    get_peak_memory()
    @printf("It took %.5f seconds\n", time()-start_time)
end

muR = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 3.0
lamda_cav = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.1
set_name = length(ARGS) >= 3 ? String(strip(ARGS[3])) : "SetII"


mus = [-3.0:1.0:4.0; 0.5:0.2:3.5]
lambdas = collect(0.0:0.025:0.2)
sets = ["SetI", "SetII", "SetIII"]

for mu in mus
    for lambda in lambdas
        for set_name in sets
            main(mu_R_val = mu, lambda_cav_val = lambda, set_name = set_name)
        end
    end
end
