include("Ising_MC_Wolff.jl")

conf = spins(0.01, -1)
Wolff_flip(conf)
println(conf.spin)
println(bond(conf, (1,2), (1,3)))
