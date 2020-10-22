include("Ising_MC.jl")
using Plots
using BenchmarkTools

"""
Calculate high energy expansion of E_site.

Parameter:
T Temperature
J Coupling constant. J = 1 is the default.

Output:
Eh::Array{Float64,3} The three entries of the output array(Eh[1], Eh[2] and Eh[3]) correspond to the result of the expansion to the order of 4,6 and 8
"""
function Eh(T, J = 1)
    Eh = []
    β = 1 / T
    E = 0
    E -= J * tanh(β * J) + 4J * tanh(β * J)^3 * (1 - tanh(β * J)^2)
    push!(Eh, E)
    E -= 2 * 6 * J * tanh(β * J)^5 * (1 - tanh(β * J)^2)
    push!(Eh, E)
    E -= 4.5 * 8 * J * tanh(β * J)^7 * (1 - tanh(β * J)^2)
    push!(Eh, E)
    return Eh
end

function EL(T, J = 1)
    EL = []
    β = 1 / T
    E = 0
    E += -2J + 2J * 4exp(-2β * J)^4
    push!(EL, E)
    E += 2J * 6 * exp(-2β * J)^6
    push!(EL, E)
    E += 2J * 8 * exp(-2β * J)^8
    push!(EL, E)
    return EL
end

function driver(dimension)
    Temp = range(0.1, 100, length = 100)
    E_site = []
    S_site = []
    B = []
    χ = []
    C = []
    for T in Temp
        mic_state = spins(T, -1, 10^6, dimension)
        mac_state = macrostate(mic_state, 10000)
        push!(E_site, mac_state.energy_site)
        push!(S_site, mac_state.polar_site)
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

function graphing()
    for i = 1:3
        Temp, E_site, S_site, B, χ, C, dimension = Result[i][1],
        Result[i][2],
        Result[i][3],
        Result[i][4],
        Result[i][5],
        Result[i][6],
        Result[i][7]

        Eh4(T) = Eh(T)[1]
        Eh6(T) = Eh(T)[2]
        Eh8(T) = Eh(T)[3]

        EL4(T) = EL(T)[1]
        EL6(T) = EL(T)[2]
        EL8(T) = EL(T)[3]

        println("It's done.")
        E_site_plt = plot(
            Temp,
            E_site,
            label = "E_site",
            xlabel = "T",
            ylabel = "E_site",
            line = (1, 1.5),
            dpi = 800,
        )
        plot!(Eh4, 4, 100, label = "High Temp expansion 4th order", linestyle = :dash)
        plot!(Eh6, 4, 100, label = "High Temp expansion 6th order", linestyle = :dash)
        plot!(Eh8, 4, 100, label = "High Temp expansion 8th order", linestyle = :dash)
        plot!(EL4, 0.01, 4, label = "Low Temp expansion 4th order", linestyle = :dash)
        plot!(EL6, 0.01, 4, label = "Low Temp expansion 6th order", linestyle = :dash)
        plot!(EL8, 0.01, 4, label = "Low Temp expansion 8th order", linestyle = :dash)
        png("E_site_$dimension")

        plot(
            Temp,
            S_site,
            label = "S_site",
            xlabel = "T",
            ylabel = "S_site",
            line = (1, 1.5),
            dpi = 800,
        )
        png("S_site_$dimension")
        plot(Temp, B, label = "B", xlabel = "T", ylabel = "B", line = (1, 1.5), dpi = 800)
        png("B_$dimension")
        plot(
            Temp,
            χ,
            label = "χ_site",
            xlabel = "T",
            ylabel = "χ",
            line = (1, 1.5),
            dpi = 800,
        )
        png("χ_$dimension")
        plot(Temp, C, label = "C", xlabel = "T", ylabel = "C", line = (1, 1.5), dpi = 800)
        png("C_$dimension")

    end
end

graphing()
