include("Ising_MC.jl")
using Plots
using BenchmarkTools


Temp = range(0.01, 10, length = 100)
E_site = []
S_site = []
B = []
χ = []
C = []
for T in Temp
    mic_state = spins(T, -1, 10^6)
    mac_state = macrostate(mic_state, 10000)
    push!(E_site, mac_state.energy_site)
    push!(S_site, mac_state.polar_site)
    push!(B, (0.5 * (mac_state.M4 / mac_state.M2^2 - 1)))
    push!(
        χ,
        (1 / mic_state.T * (mac_state.M2 - mac_state.abs_polar^2)) / mic_state.dimension^2,
    )
    push!(
        C,
        1 / mic_state.T *
        (mac_state.E2_site - mac_state.energy_site^2) *
        mic_state.dimension^2,
    )
end
println("It's done.")
plot(
    Temp,
    E_site,
    label = "E_site",
    xlabel = "T",
    ylabel = "E_site, S_site, B, χ, C",
    line = (1, 1.5)
)
plot!(Temp, S_site, label = "S_site", line = (1, 1.5))
plot!(Temp, B, label = "B", line = (1, 1.5))
plot!(Temp, χ, label = "χ_site", line = (1, 1.5))
plot!(Temp, C, label = "C", line = (1, 1.5))
