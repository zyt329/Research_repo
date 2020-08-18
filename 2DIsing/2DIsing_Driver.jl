module IsingDriver

using 2DIsing
using Plots
using BenchmarkTools

import Ising_MC.spins, Ising_MC.macrostate

Temp = range(0.01, 10, length = 100)
E_site = []
for T in Temp
    push!(E_site, macrostate(spins(T,-1, 10^6),10000).energy_site)
end
plot(Temp, E_site)

end
