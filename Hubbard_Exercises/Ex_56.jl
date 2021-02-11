using LinearAlgebra
using Plots

"""Dispersion relation for up spins"""
E_up(n_down, k, U; t = 1) = -2 * t * cos(k) + U * n_down
"""Dispersion relation for down spins"""
E_down(n_up, k, U; t = 1) = -2 * t * cos(k) + U * n_up

"""finding the minimum energy amoung different magnetizations."""
function meanfield(ρ, U; N::Int = 10^3, t = 1, E_up = E_up, E_down = E_down)
    @assert(ρ <= 1, "Density should be smaller than 1.")
    N_tot = Int(ceil(ρ * N))
    #Generate k points
    kpoints = range(-π, π, length = N+1)
    E_tots = zeros(N_tot)

    for N_up = 1:N_tot
        N_down = N_tot - N_up
        for k in kpoints
            if k >= ((-N_up / 2 + 1) * 2π / N) && k <= (N_up / 2 * 2π / N)
                E_tots[N_up] += E_up(N_down / N, k, U; t = t)
            end
        end
        if isodd(N_up)
            E_tots[N_up] += E_up(N_down / N, (N_up / 2 +1/2)* 2π / N, U; t = t)
        end
        for k in kpoints
            if k >= ((-N_down / 2 + 1) * 2π / N) && k <= (N_down / 2 * 2π / N)
                E_tots[N_up] += E_down(N_up / N, k, U; t = t)
            end
        end
        if isodd(N_up)
            E_tots[N_up] += E_down(N_up / N, (N_down / 2 +1/2)* 2π / N, U; t = t)
        end
        E_tots[N_up] = E_tots[N_up] / N - U * N_up / N * N_down / N
    end


    N_up_best = argmin(E_tots)
    return (
        E_tots[N_up_best],
        N_up_best,
        (2 * N_up_best - N_tot) / N_tot,
        plot(1:N_tot, E_tots),
    )
end

meanfield(3/4, 5.163)[4]
