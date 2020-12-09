using JLD
using Plots

Result = load("Result_Thu_03_Dec_2020_01_35_45.jld")["Result"]
#=
Result = load("Result.jld")["Result"]
push!(Result,load("Result4_Mon_30_Nov_2020_21_53_51.jld")["Result"])
=#

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
    E -= 2 * J * tanh(β * J) + 4J * tanh(β * J)^3 * (1 - tanh(β * J)^2)
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

function graphing()
    for i = 1:4
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
            Temp[1:Int(floor(0.8 * end))],
            E_site[1:Int(floor(0.8 * end))],
            label = "E_site",
            xlabel = "T",
            ylabel = "E_site",
            line = (1, 1.5),
            dpi = 800,
            legend = :topleft,
        )
        plot!(Eh4, 2, 5, label = "High Temp expansion 4th order", linestyle = :dash)
        plot!(Eh6, 2, 5, label = "High Temp expansion 6th order", linestyle = :dash)
        plot!(Eh8, 2, 5, label = "High Temp expansion 8th order", linestyle = :dash)
        plot!(EL4, 0.01, 4, label = "Low Temp expansion 4th order", linestyle = :dash)
        plot!(EL6, 0.01, 4, label = "Low Temp expansion 6th order", linestyle = :dash)
        plot!(EL8, 0.01, 4, label = "Low Temp expansion 8th order", linestyle = :dash)
        gui(E_site_plt)
        savefig("E_site_$dimension.svg")

        plot(
            Temp[1:Int(floor(0.8 * end))],
            S_site[1:Int(floor(0.8 * end))],
            label = "S_site",
            xlabel = "T",
            ylabel = "S_site",
            line = (1, 1.5),
            dpi = 800,
        )
        savefig("S_site_$dimension.svg")

        plot(
            Temp[1:Int(floor(0.8 * end))],
            B[1:Int(floor(0.8 * end))],
            label = "B",
            xlabel = "T",
            ylabel = "B",
            line = (1, 1.5),
            dpi = 800,
        )
        savefig("B_$dimension.svg")

        plot(
            Temp[1:Int(floor(0.8 * end))],
            χ[1:Int(floor(0.8 * end))],
            label = "Chi_site",
            xlabel = "T",
            ylabel = "Chi",
            line = (1, 1.5),
            dpi = 800,
        )
        plot!([peaks[1][i][2], peaks[1][i][2]], [peaks[1][i][1], 0])
        savefig("χ_$dimension.svg")
        plot(
            Temp[1:Int(floor(0.8 * end))],
            C[1:Int(floor(0.8 * end))],
            label = "C",
            xlabel = "T",
            ylabel = "C",
            line = (1, 1.5),
            dpi = 800,
        )
        plot!([peaks[2][i][2], peaks[2][i][2]], [peaks[2][i][1], 0])
        savefig("C_$dimension.svg")
    end
end

function graphS(Result)
    for i = 1:4
        plot(
            Result[i][1][1:Int(floor(0.8 * end))],
            Result[i][3][1:Int(floor(0.8 * end))],
            label = "S_site",
            xlabel = "T",
            ylabel = "S_site",
            line = (1, 1.5),
            dpi = 800,
        )
        plot!(
            Result[i][1][1:Int(floor(0.8 * end))],
            abs.(Result[i][3][1:Int(floor(0.8 * end))]),
            label = "S_site_abs",
            line = (1, 1.5),
            linestyle = :dash,
        )
        savefig("S_site_all$(Result[i][7])")
    end
end
function Bgraphing()
    plot()
    for i = 1:4
        Temp, E_site, S_site, B, χ, C, dimension = Result[i][1],
        Result[i][2],
        Result[i][3],
        Result[i][4],
        Result[i][5],
        Result[i][6],
        Result[i][7]

        plot!(
            Temp[Int(floor(0.3 * end)):Int(floor(0.4 * end))],
            B[Int(floor(0.3 * end)):Int(floor(0.4 * end))],
            label = "B_$dimension",
            xlabel = "T",
            ylabel = "B",
            line = (1, 1.5),
            dpi = 800,
            legend = :topleft,
        )
    end
    plot!([2.26, 2.26], [0, 0.5], linestyle = :dash, label = "")
    savefig("B_crossing.svg")
end
Bgraphing()
#graphing()
#graphS(Result)
