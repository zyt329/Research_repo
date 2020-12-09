# Analyzation of critical exponents
using Plots
using BenchmarkTools
using JLD
using Dates
using LsqFit

peaks = load("peaks_Wed_02_Dec_2020_17_57_11.jld")["peaks"]
Result = load("Result_Thu_03_Dec_2020_01_35_45.jld")["Result"]
#=
Result = load("Result.jld")["Result"]
push!(Result, load("Result4_Wed_02_Dec_2020_16_03_35.jld")["Result"])
=#

"""
Plot Tc(L) versus L^-1 using data for χ.

"""
function ex4(peaks)
    L = [8, 12, 16, 24] .^ (-1)
    Tc = []
    for i = 1:4
        push!(Tc, peaks[1][i][2])
    end
    plot(
        L,
        Tc,
        seriestype = :scatter,
        markerstrokewidth = 0.3,
        markersize = 2,
        title = "TC(L) vs L^-1 for Chi",
        xlabel = "L^-1",
        ylabel = "TC",
        label = "Tc",
        legend = :topleft,
        dpi = 800,
    )


    model(L, P) = P[1] * L .+ P[2]
    fit = curve_fit(model, L, Tc, [5, 2.5])
    println(fit.param)

    y(x) = 2.105x + 2.271
    plot!(y, 0, 0.15, label = "Fitted line")
    savefig("TC(L) vs L^-1 for E_site_16.png.svg")

end
#The critical temperature is 2.2719 according to the simulation results.
#println(peaks[1][4][2])
ex4(peaks)

function ex5_1(peaks)
    logL = log.([8, 12, 16, 24])
    χ = []
    for i = 1:4
        push!(χ, peaks[1][i][1])
    end
    logχ = log.(χ)
    plot(
        logL,
        logχ,
        seriestype = :scatter,
        markerstrokewidth = 0.3,
        markersize = 2,
        dpi = 800,
        xlabel = "logL",
        ylabel = "logChi(Tc(L))",
        label = "logChi(Tc(L))",
        title = "logChi(Tc(L)) vs logL",
        legend = :topleft
    )

    y(x) = 1.8499094765645645x -3.326866209162157
    plot!(y, 2, 3.25, label = "Fitted line")
    savefig("EX5.svg")

    logχmodel(logL, P) = P[1] * logL .+ P[2]
    fit = curve_fit(logχmodel, logL, logχ, [1.0, 1.0])
    println(fit.param)
end
ex5_1(peaks)
# γ/ν is 1.8499 according to my simulation results.1.8499094765645645, -3.326866209162157

function ex5_2(Result, peaks)
    plot(dpi = 800)
    for i = 1:4
        L = Result[i][7]
        χ = Result[i][5][1:Int(floor(0.8 * end))]
        TcL = peaks[1][i][2]
        Temp = Result[i][1][1:Int(floor(0.8 * end))]

        yaxis = L^(-1.8499) .* χ
        xaxis = L .* (Temp .- TcL)

        plot!(
            xaxis,
            yaxis,
            seriestype = :scatter,
            markerstrokewidth = 0.3,
            markersize = 2,
            title = "Collpase of scaled Chi",
            xlabel = "L^1/nu (T - Tc(L))",
            ylabel = "L^(gamma/nu)Chi",
            label = "Scaled Chi(L=$L)",
        )

    end
    savefig("EX5_2.svg")
end
ex5_2(Result, peaks)

function ex6(Result, peaks)
    plot(dpi = 800)
    for i = 1:4
        L = Result[i][7]
        S = Result[i][3][1:Int(floor(0.8 * end))]
        Temp = Result[i][1][1:Int(floor(0.8 * end))]

        yaxis = L^(0.25) * S * L^(-2)
        xaxis = Temp
        plot!(
            xaxis,
            yaxis,
            title = "Intersection of M",
            xlabel = "T",
            ylabel = "L^(beta/nu)M",
            label = "L^(beta/nu)M(L=$L)",
        )

    end
    savefig("Intersection of M.svg")
    #The intersection is at Tc'=2.32
end
ex6(Result, peaks)

function ex6_2(Result, peaks)
    plot(dpi = 800)
    for i = 1:4
        L = Result[i][7]
        S = Result[i][3][1:Int(floor(0.8 * end))]
        Temp = Result[i][1][1:Int(floor(0.8 * end))]

        yaxis = L^(0.25) * S * L^(-2)
        xaxis = L * (Temp .- 2.32)
        plot!(
            xaxis,
            yaxis,
            title = "Collapse of L^(beta/nu)M",
            xlabel = "L^1/nu (T-Tc')",
            ylabel = "L^(beta/nu)M",
            label = "L^(beta/nu)M(L=$L)",
            seriestype = :scatter,
            markerstrokewidth = 0.3,
            markersize = 2,
        )

    end
    savefig("Collapse of Scaled M.svg")
end
ex6_2(Result, peaks)

function ex_7_1(Result, peaks)
    L = [8, 12, 16, 24] .^ (-1)
    Tc = []
    C = []
    for i = 1:4
        push!(C, peaks[2][i][1])
        push!(Tc, peaks[2][i][2])
    end
    plot(
        L,
        Tc,
        seriestype = :scatter,
        markerstrokewidth = 0.3,
        markersize = 2,
        title = "TC(L) vs L^-1 for C",
        xlabel = "L^-1",
        ylabel = "TC",
        label = "Tc",
        legend = :topleft,
        dpi = 800,
    )
    y(x) = 1.313x + 2.2508
    plot!(y, 0, 0.15, label = "Fitted line")

    savefig("ex7_1.svg")

    Tcmodel(L, P) = P[1] * L .+ P[2]
    fit = curve_fit(Tcmodel, L, Tc, [1.0, 1.0])
    println(fit.param)
end
ex_7_1(Result, peaks)


function ex_7_2(Result, peaks)
    L = [8, 12, 16, 24] .^ (-1)
    Tc = []
    C = []
    for i = 1:4
        push!(C, peaks[2][i][1])
        push!(Tc, peaks[2][i][2])
    end
    logL = log.([8, 12, 16, 24])
    plot(
        logL,
        C,
        seriestype = :scatter,
        markerstrokewidth = 0.3,
        markersize = 2,
        title = "Cmax(L) vs Log(L)",
        xlabel = "Log(L)",
        ylabel = "Cmax",
        label = "Cmax",
        legend = :topleft,
        dpi = 800,
    )
    Cmaxmodel(logL, P) = P[1] * logL .+ P[2]
    fit = curve_fit(Cmaxmodel, logL, C, [1.0, 1.0])
    println(fit.param)
    y(x) = 0.924x + 0.9668
    plot!(y, 2, 3.25, label = "Fitted line")
    savefig("ex7_2.svg")
end
ex_7_2(Result, peaks)
