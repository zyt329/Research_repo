include("Ising_MC.jl")
#module Ising_Driver
using Plots
using BenchmarkTools

#import Main.Ising_MC.spins, Main.Ising_MC.macrostate

Temp = range(0.01, 10, length = 100)
E_site = []
S_site = []
for T in Temp
    mac_state = macrostate(spins(T,-1, 10^5),10000)
    push!(E_site, mac_state.energy_site)
    push!(S_site, mac_state.polar_site)
end
#println("It's done.")
plot(Temp, E_site)
plot!(Temp, S_site)
