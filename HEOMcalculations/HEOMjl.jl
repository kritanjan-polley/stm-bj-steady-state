using HierarchicalEOM
using Printf

struct PROCESS_MEMORY_COUNTERS
    cb::UInt32
    PageFaultCount::UInt32
    PeakWorkingSetSize::UInt
    WorkingSetSize::UInt
    QuotaPeakPagedPoolUsage::UInt
    QuotaPagedPoolUsage::UInt
    QuotaPeakNonPagedPoolUsage::UInt
    QuotaNonPagedPoolUsage::UInt
    PagefileUsage::UInt
    PeakPagefileUsage::UInt
end

function get_peak_memory_bytes()::Int64
    if Sys.islinux() || Sys.isapple()
        rusage = zeros(Int64, 18)
        ret = ccall(:getrusage, Int32, (Int32, Ptr{Cvoid}), 0, rusage)

        if ret == 0
            # rusage[1-2] = utime, rusage[3-4] = stime, rusage[5] = maxrss
            multiplier = Sys.islinux() ? 1024 : 1
            return rusage[5] * multiplier
        end

    elseif Sys.iswindows()
        hProcess = ccall(:GetCurrentProcess, Ptr{Cvoid}, ())

        mem_counters = Ref(PROCESS_MEMORY_COUNTERS(0,0,0,0,0,0,0,0,0,0))
        cb = sizeof(PROCESS_MEMORY_COUNTERS)

        ret = ccall((:GetProcessMemoryInfo, "psapi"), Int32,
                    (Ptr{Cvoid}, Ptr{PROCESS_MEMORY_COUNTERS}, UInt32),
                    hProcess, mem_counters, cb)

        if ret != 0
            return Int64(mem_counters[].PeakWorkingSetSize)
        end
    end

    return -1
end

start_time = time()
println("Initial Peak: $(get_peak_memory_bytes() / (1024^3)) GB")
println("Starting simulation")

set_name = length(ARGS) >= 3 ? strip(ARGS[3]) : "SetI"
epsilon1 = set_name == "SetI" ? -1.0 : -0.5
epsilon2 = -epsilon1
jval = 0.0
factor_omega = set_name == "SetII" ? 0.4 : 1.0
omega = factor_omega * abs(epsilon2 - epsilon1)
Gamma = 1.0
lamda = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 0.1
W = 10.0
muL = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) * -1.0 : -1.0
muR = -muL
kBT = 1.0
Np = 8
tmax = 20.0
time_step = 0.1
lamda_dir = @sprintf("lamda_%.4f", lamda)
mu_dir = @sprintf("mu_%.1f", isapprox(muL, 0.0) ? 0.0 : muL)

if length(ARGS) >= 2
    filestr = (
        "$(lamda_dir)_$(mu_dir)"
        )
else
    # filestr = "test$(Np)"
    filestr = "test2"
end

println("Parameters for set: $(set_name)")
println("File prefix: $(filestr)")
@printf("Delta (diabatic energy bias) = %.4f\n", abs(epsilon1 - epsilon2))
@printf("Diabatic coupling = %.4f\n", jval)
@printf("Oscillator frequency = %.4f\n", omega)
@printf("Gamma = %.4f\n", Gamma)
@printf("Lamda = %.4f\n", lamda)
@printf("Fermi energy (left) = %.4f\n", muL)
@printf("Fermi energy (right) = %.4f\n", muR)
@printf("kBT = %.4f\n", kBT)
@printf("Band width = %.4f\n", W)
@printf("Maximum propagation time = %.4f\n", tmax)
@printf("Time Step = %.4f\n", time_step)

a_boson = destroy(Np + 1)
adag_boson = create(Np + 1)
eye_boson = qeye(Np + 1)
d_ferm = destroy(2)
ddag_ferm = create(2)
d_eye = qeye(2)

d1_coup = tensor(d_ferm, d_eye, eye_boson)
d2_coup = tensor(d_eye, d_ferm, eye_boson)

hs = (
    epsilon1 * tensor(ddag_ferm * d_ferm, d_eye, eye_boson)
    + epsilon2 * tensor(d_eye, ddag_ferm * d_ferm, eye_boson)
    + jval * (
        tensor(ddag_ferm, d_ferm, eye_boson)
        + tensor(d_ferm, ddag_ferm, eye_boson)
        )
    )
hp = omega * tensor(d_eye, d_eye, (adag_boson * a_boson
                    + qeye(Np + 1) * 0.5))

hsp = (lamda * (
    tensor(ddag_ferm, d_ferm, eye_boson) +
    tensor(d_ferm, ddag_ferm, eye_boson)
    ) * tensor(d_eye, d_eye, a_boson + adag_boson)
    )

h = hs + hp + hsp
psi0 = tensor(basis(2, 0), basis(2, 0), basis(Np + 1, 0))

# println(h)

N = 3
nbath = 2
bath_L1 = Fermion_Lorentz_Pade(d1_coup, Gamma, muL, W, kBT, N)
bath_L2 = Fermion_Lorentz_Pade(d2_coup, Gamma, muL, W, kBT, N)
bath_R1 = Fermion_Lorentz_Pade(d1_coup, Gamma, muR, W, kBT, N)
bath_R2 = Fermion_Lorentz_Pade(d2_coup, Gamma, muR, W, kBT, N)
if nbath == 4
    baths = [bath_L1, bath_L2, bath_R1, bath_R2]
else
    baths = [bath_L1, bath_R2]
end

tier = 4 # 4
@assert length(baths) == 4 || length(baths) == 2
println("The system modes have $(length(baths)) connected baths")
println("Using $(Np+1) levels for the vibrational mode")
println("Using $(N) exponential terms for spectral density")
println("Using $(tier) levels of hierarchy")
M = M_Fermion(h, tier, baths; threshold = 1e-7)
println(M)
time_arr = range(0.0, tmax, step = time_step)
ados_evolution = HEOMsolve(M, psi0, time_arr).ados

println("Done with HEOM evolution")

# operators
println("Getting expectation values")
time_here = time()
photon_num_op = tensor(d_eye, d_eye, adag_boson * a_boson)
pop1_op = tensor(ddag_ferm * d_ferm, d_eye, eye_boson)
pop2_op = tensor(d_eye, ddag_ferm * d_ferm, eye_boson)

function Ic(ados, M::M_Fermion, bathIdx::Int)
    HDict = M.hierarchy

    idx_list = HDict.lvl2idx[1]
    I = 0.0im
    for idx in idx_list
        rho1 = ados[idx]

        nvec = HDict.idx2nvec[idx]
        for (alpha, k, _) in getIndexEnsemble(nvec, HDict.bathPtr)
            if alpha == bathIdx
                exponent = M.bath[alpha][k]
                if exponent.types == "fA"
                    I += tr(exponent.op' * rho1)
                elseif exponent.types == "fE"
                    I -= tr(exponent.op' * rho1)
                end
                break
            end
        end
    end

    return real(1im * I)
end

# expectation values
photon_number = expect(photon_num_op, ados_evolution)
pop1_expect = expect(pop1_op, ados_evolution)
pop2_expect = expect(pop2_op, ados_evolution)

N_b = Np + 1
boson_projectors = [tensor(d_eye, d_eye, basis(N_b, n) * basis(N_b, n)') for n in 0:Np]
boson_populations = zeros(length(time_arr), N_b)

for (ti, ados) in enumerate(ados_evolution)
    rho = ados[1]
    for (n, Pn) in enumerate(boson_projectors)
        boson_populations[ti, n] = real(tr(Pn * rho))
    end
end

Ie_L = Float64[]
Ie_R = Float64[]
for ados in ados_evolution
    if length(baths) == 2
        left1 = Ic(ados, M, 1)
        right2 = Ic(ados, M, 2)
        push!(Ie_L, left1)
        push!(Ie_R, right2)
    elseif length(baths) == 4
        left1 = Ic(ados, M, 1)
        left2 = Ic(ados, M, 2)
        right1 = Ic(ados, M, 3)
        right2 = Ic(ados, M, 4)
        push!(Ie_L, left1 + left2)
        push!(Ie_R, right1 + right2)
    else
        # assert should prevent this
        println("Something is wrong with the bath modes")
        exit()
    end
end


@printf("%.5f seconds for expectation values\n", time()-time_here)
println("Peak: $(get_peak_memory_bytes() / (1024^3)) GB")

photon_number_file = joinpath(set_name, lamda_dir, mu_dir, "photon.HEOM.txt")
mkpath(dirname(photon_number_file))
f = open(photon_number_file, "w")
for (i, val) in enumerate(time_arr)
    @printf(f, "%10.6f %10.6f\n", val, photon_number[i])
end
close(f)

population_file = joinpath(set_name, lamda_dir, mu_dir, "population.HEOM.txt")
mkpath(dirname(population_file))
f = open(population_file, "w")
for (i, val) in enumerate(time_arr)
    @printf(f, "%10.6f %10.6f %10.6f\n",
            val, pop1_expect[i], pop2_expect[i])
end
close(f)

current_file = joinpath(set_name, lamda_dir, mu_dir, "current.HEOM.txt")
mkpath(dirname(current_file))
f = open(current_file, "w")
for (i, val) in enumerate(time_arr)
    @printf(f, "%12.6f  %12.6f  %12.6f\n", val, Ie_L[i], Ie_R[i])
end
close(f)

boson_file = joinpath(set_name, lamda_dir, mu_dir, "boson_diag.HEOM.txt")
mkpath(dirname(boson_file))
open(boson_file, "w") do f
    for (i, t) in enumerate(time_arr)
        data_str = join([@sprintf("%10.8f", x) for x in boson_populations[i, :]], " ")
        println(f, @sprintf("%10.6f ", t) * data_str)
    end
end
close(f)

println("Memory at the end: $(get_peak_memory_bytes() / (1024^3)) GB")
@printf("It took %.5f seconds\n", time()-start_time)
