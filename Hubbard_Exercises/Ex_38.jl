using Plots

function Gf1D(τ; T = 10, n::Int64 = 0, l::Int64 = 1, t = 1, E0 = 2, N = 10^3)
    E(k) = E0 - 2t * cos(k)
    f(k) = 1 / (exp(E(k) / T + 1))
    kpoints = range(-π, π, length = N)
    function G(t)
        G = 0
        for k in kpoints
            G += 1 / N * exp( im * k * (n - l)) * (1 - f(k)) * exp(-E(k) * t)
        end
        G
    end
    #=function ReG(t)
        ReG = 0
        for k in kpoints
            ReG += 1 / N * sin(k * (n - l)) * (1 - f(k)) * exp(-E(k) * t)
        end
        ReG
    end=#
    return G(τ)
end

function Gf2D(τ; T = 10, n = [0, 0], l = [1, 0], t = 1, E0 = 2, N = 10^2)
    E(k) = -2t * (cos(k[1]) + cos(k[2]))
    f(k) = 1 / (exp(E(k) / T + 1))
    kpoints = []
    for kx in range(-π, π, length = N)
        for ky in range(-π, π, length = N)
            push!(kpoints, [kx, ky])
        end
    end
    function G(t)
        G = 0
        for k in kpoints
            G += 1 / N^2 * exp( im * dot(k, (n - l))) * (1 - f(k)) * exp(-E(k) * t)
        end
        return G
    end
    #=
    function ReG(t)
        ReG = 0
        for k in kpoints
            ReG += 1 / N * sin(dot(k, (n - l))) * (1 - f(k)) * exp(-E(k) * t)
        end
        return ReG
    end
    =#
    return G(τ)
end

Gf1DIm(τ) = imag(Gf1D(τ, T = 0.1))
Gf1DRe(τ) = real(Gf1D(τ, T = 0.1))
p1 = plot(Gf1DIm, -1, 1, label = "Imaginary Part", title = "1D")
plot!(Gf1DRe, -1, 1, label = "Real Part")

Gf2DIm(τ) = imag(Gf2D(τ, T = 0.1))
Gf2DRe(τ) = real(Gf2D(τ, T = 0.1))
p2 = plot(Gf2DIm, -1, 1, label = "Imaginary Part", title = "2D")
plot!(Gf2DRe, -1, 1, label = "Real Part")

plot(p1, p2, layout=(2,1))
