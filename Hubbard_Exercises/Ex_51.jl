using Plots
using LinearAlgebra

function oneuponedown(u, t, T, μ = 0)
    H = vcat(
        hcat(-1 / 2 * u - 2 * μ, -2 * t, -t, t),
        hcat(-2t, -1 / 2 * u - 2μ, -t, -t),
        hcat(-t, -t, 1 / 2 * u - 2μ, -2t),
        hcat(-t, -t, -2t, 1 / 2 * u - 2μ),
    )
    M2 = [
        1 0 0 0
        0 1 0 0
        0 0 0 0
        0 0 0 0
    ]
    eigH = eigen(H)
    m2 = 0
    ρ = zeros(4, 4)
    for k = 1:length(eigH.values)
        ρ[k, k] = exp(-eigH.values[k] / T)
    end
    ρ = ρ / tr(ρ)
    include("Ex_51debug.jl")
    tr(ρ * eigH.vectors' * M2 * eigH.vectors)
end

M2u(u) = oneuponedown(u, 1, 1)
p1 = plot(M2u, 0, 100, xlabel = "U", label = "U", color = "red")
M2t(t) = oneuponedown(1, t, 1)
p2 = plot(M2t, 0, 1, xlabel = "t", label = "t", color= "blue")
M3T(T) = oneuponedown(1, 1, T)
p3 = plot(M3T, 0, 100, xlabel = "T", label = "T",color = "black")
plot(p1, p2, p3, layout = (3, 1))
