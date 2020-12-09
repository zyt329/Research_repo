include("Ising_MC.jl")
using Plots
using BenchmarkTools
using JLD
using Dates

function driver(dimension)
    Temp = vcat(range(0.1, 5, length = 400),range(5.01,100, length = 100))
    E_site = []
    S_site = []
    B = []
    χ = []
    C = []
    for T in Temp
        mic_state = spins(T, -1, Int(floor(10^6
        *(dimension/8)^2)), dimension)
        mac_state = macrostate(mic_state, 10000)
        push!(E_site, mac_state.energy_site)
        #I used absolute value of S in the place of S
        push!(S_site, mac_state.abs_polar)
        push!(B, (0.5 * (mac_state.M4 / mac_state.M2^2 - 1)))
        push!(
            χ,
            (1 / mic_state.T * (mac_state.M2 - mac_state.abs_polar^2)) /
            mic_state.dimension^2,
        )
        push!(
            C,
            1 / mic_state.T *
            (mac_state.E2_site - mac_state.energy_site^2) *
            mic_state.dimension^2,
        )
    end
    return (Temp, E_site, S_site, B, χ, C, dimension)
end

Result = []
push!(Result, driver(8))
push!(Result, driver(12))
push!(Result, driver(16))
push!(Result, driver(24))


my_time = Dates.now()
save("Result_$(Dates.format(my_time, "e_dd_u_yyyy_HH_MM_SS")).jld", "Result", Result)
println("Simulation is completed")
#=
Result4 = driver(24)
my_time = Dates.now()
save("Result4_$(Dates.format(my_time, "e_dd_u_yyyy_HH_MM_SS")).jld", "Result", Result4)
=#
