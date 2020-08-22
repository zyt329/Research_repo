
"""
Represent a state.

fields:  T, J, steps, dimension, spin
"""
mutable struct spins
    T::Any
    J::Any
    steps::Int64
    dimension::Int64
    spin::Array{Int64,2}
    energy::Any
    polar::Any
    function spins(
        T,
        J,
        steps::Int64 = 10^5,
        dimension::Int64 = 8,
        spin::Array{Int64,2} = ones(Int64, dimension, dimension),
    )
        S = spin
        N = dimension
        E0 = 0
        for i = 1:N
            for j = 1:N
                E0 =
                    E0 +
                    0.5 *
                    J *
                    S[i, j] *
                    (
                        S[mod1(i - 1, N), j] +
                        S[mod1(i + 1, N), j] +
                        S[i, mod1(j - 1, N)] +
                        S[i, mod1(j + 1, N)]
                    )
            end
        end
        polar = sum(spin)
        new(T, J, steps, dimension, spin, E0, polar)
    end
end

"""
Calculate the energy difference after flipping one spin.

Parameter:
conf::spins : Input the state
index::Tuple : An Tuple of index of the flipped spin.
"""
function Energy_Diff(conf::spins, index::Tuple)
    N = conf.dimension
    S = conf.spin
    m = index[1]
    n = index[2]
    @assert(
        m <= N && m >= 1 && n <= N && n >= 1,
        "The Index of flipped spin should be within range"
    )
    ΔE =
        -conf.J *
        2 *
        S[m, n] *
        (
            S[mod1(m - 1, N), n] +
            S[mod1(m + 1, N), n] +
            S[m, mod1(n - 1, N)] +
            S[m, mod1(n + 1, N)]
        )
    return ΔE
end

"""
Choose a spin to flip

Parameter:
conf::spins : the system

Return:
A tuple of index of the spin to flip
"""
function flip_index(conf::spins)
    (rand(1:conf.dimension), rand(1:conf.dimension))
end

"""
Function to flip a spin

Parameter:
conf::spins : current configuration
index::Tuple : indexs of the flipped spin

Return:
new_config::spins : the new configuration
"""
function flip(conf::spins, index::Tuple)
    conf.energy += Energy_Diff(conf, index)
    conf.polar += -2 * conf.spin[index[1], index[2]]
    conf.spin[index[1], index[2]] = -1 * conf.spin[index[1], index[2]]
    return conf
end

"""
Calculating the energy of the current state

Parameter:
conf::spins : the current state

Return:
Energy::Float64 : energy of the current state
"""
function Energy(conf::spins)
    S = conf.spin
    N = conf.dimension
    E0 = 0
    for i = 1:N
        for j = 1:N
            E0 =
                E0 +
                0.5 *
                conf.J *
                S[i, j] *
                (
                    S[mod1(i - 1, N), j] +
                    S[mod1(i + 1, N), j] +
                    S[i, mod1(j - 1, N)] +
                    S[i, mod1(j + 1, N)]
                )
        end
    end
    return E0
end

"""
Deciding whether to flip a chosen spin

Parameter:
conf::spins : The current state.
index::Tuple : A tuple of spin index to flip.

Return:
true of false. true means to flip false means not to
"""
function evol_cond(conf::spins, index::Tuple)
    #index = flip_index(conf.dimension)
    ΔE = Energy_Diff(conf, index)

    p = exp(-ΔE / conf.T) / (1 + exp(-ΔE / conf.T))
    r = rand()
    if r < p
        return true
    else
        return false
    end
end

"""
Represent a macrostate.

Fields:
energy_site::Float64 : Energy per site for the system at a certain temperature.
polar_site::Float64 : Polarization per site for the system at a certain temperature.

"""
struct macrostate
    energy_site::Float64
    E2_site::Float64
    polar_site::Float64
    abs_polar::Float64
    M2::Float64
    M4::Float64
    """
    One method of the constructor's for the macrostate type.

    Parameter:
    Conf::spins : The lattice system
    cutoff::Int64 : Determines to how many steps would the system arrive at equilibrium(We would only calculate quantities after the cutoff step)
    """
    function macrostate(conf::spins, cutoff::Int64)
        energy_site = 0
        E2_site = 0
        polar_site = 0
        abs_polar = 0
        M2 = 0
        M4 = 0
        for i = 1:conf.steps
            index = flip_index(conf)
            if evol_cond(conf, index)
                flip(conf, index)
            end
            i < cutoff + 1 && continue

            energy_site += conf.energy / conf.dimension^2
            E2_site += (conf.energy / conf.dimension^2)^2
            polar_site += conf.polar / conf.dimension^2
            abs_polar += abs(conf.polar)
            M2 += conf.polar^2
            M4 += conf.polar^4
        end
        energy_site = energy_site / (conf.steps - cutoff)
        E2_site = E2_site / (conf.steps - cutoff)
        polar_site = polar_site / (conf.steps - cutoff)
        abs_polar = abs_polar / (conf.steps - cutoff)
        M2 = M2 / (conf.steps - cutoff)
        M4 = M4 / (conf.steps - cutoff)

        new(energy_site, E2_site, polar_site, abs_polar, M2, M4)
    end
end
