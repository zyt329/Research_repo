include("Ising_MC.jl")
using Plots
using BenchmarkTools
using JLD
using Dates

function driver(dimension)
    #Temp = vcat(range(0.1, 5, length = 400),range(5.01,100, length = 100))
    Temp = (2.1, 2.2, 2.3, 2.4)
    E_site = Float64[]
    S_site = Float64[]
    B = Float64[]
    χ = Float64[]
    C = Float64[]
    mac_states = macrostate{Array{Float64,1}}[]
    for T in Temp
        #mic_state = spins(T, -1, Int(floor(10^6*(dimension/8)^2)), dimension)
        mic_state = spins(T, -1, Int(floor(10^5)), dimension)
        mac_state = macrostate(mic_state, 10000)
        averages = avgs(mac_state)
        push!(mac_states, mac_state)
        push!(E_site, averages.energy_site)
        #I used absolute value of S in the place of S
        push!(S_site, averages.abs_polar)
        #push!(B, (0.5 * (mac_state.M4 / mac_state.M2^2 - 1)))
        push!(B, (1 - averages.M4 /(3* averages.M2^2)))
        push!(
            χ,
            (1 / mic_state.T^2 * (averages.M2 - averages.abs_polar^2)) /
            mic_state.dimension^2,
        )
        push!(
            C,
            1 / mic_state.T^2 *
            (averages.E2_site - averages.energy_site^2) *
            mic_state.dimension^2,
        )
    end
    return (Temp, E_site, S_site, B, χ, C, dimension, mac_states)
end

#Result = []
#push!(Result, driver(8))
#push!(Result, driver(12))
#push!(Result, driver(16))
#push!(Result, driver(24))

#=
my_time = Dates.now()
save("Result_$(Dates.format(my_time, "e_dd_u_yyyy_HH_MM_SS")).jld", "Result", Result)
println("Simulation is completed")
=#
#=
Result4 = driver(24)
my_time = Dates.now()
save("Result4_$(Dates.format(my_time, "e_dd_u_yyyy_HH_MM_SS")).jld", "Result", Result4)
=#
