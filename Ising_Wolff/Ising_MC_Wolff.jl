using Statistics

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
Function to decide whether to form a bond.

Input:
conf::spins : current state of the lattice
x::Tuple : A Tuple of index of spin that's inside the Wolff's cluster(flipped spin in the Wolff's algorithm)
y::Tuple : A Tuple of index of spin that represents the spin to be flipped

output:
true for a bond, false for no bond
"""
function bond(conf::spins, x, y)
    T = conf.T
    J = conf.J
    p = 1 - exp(min(0, - J / T * conf.spin[x[1], x[2]] * 2 * conf.spin[y[1], y[2]]))
    r = rand()
    if r < p
        return true
    else
        return false
    end

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
function to construct the Wolff cluster and flip all of the spins in the cluster.

Input:
conf::spins : configuration of the system.


Output:


"""
function Wolff_flip(conf::spins)
    """
    Search function to be performed recursively to complete a depth-first search of the entire lattice for the Wolff Cluster. Meanwhile flip the spins that's been reached.

    Input :
    conf::spins : configuration of the current state.
    current_index::Tuple : index of the current position.

    """
    discovered = []
    function search_flip(conf::spins, current_index::Tuple)
        N = conf.dimension
        conf = flip(conf, current_index)
        push!(discovered, current_index)
        neighbours = (
            (mod1(current_index[1] + 1, N), current_index[2]),
            (current_index[1], mod1(current_index[2] + 1, N)),
            (mod1(current_index[1] - 1, N), current_index[2]),
            (current_index[1], mod1(current_index[2] - 1, N)),
        )
        for neighb in neighbours
            if (neighb ∉ discovered) && bond(conf, current_index, neighb)
                search_flip(conf, neighb)
            end
        end
    end
    search_flip(conf, flip_index(conf))
    #println(discovered)
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


#=
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
=#

"""
An alternate to represent a macrostate. Stores all simulated states in order to do detailed analysis like error analysis or so.

Fields:
energy_site::Float64 : Energy per site for the system at a certain temperature.
polar_site::Float64 : Polarization per site for the system at a certain temperature.

"""
struct macrostate
    energy_site::Array
    E2_site::Array
    polar_site::Array
    abs_polar::Array
    M2::Array
    M4::Array
    """
    One method of the constructor's for the macrostate type.

    Parameter:
    Conf::spins : The lattice system
    cutoff::Int64 : Determines to how many steps would the system arrive at equilibrium(We would only calculate quantities after the cutoff step)
    """
    function macrostate(conf::spins, cutoff::Int64)
        energy_site = []
        E2_site = []
        polar_site = []
        abs_polar = []
        M2 = []
        M4 = []
        for i = 1:conf.steps
            Wolff_flip(conf)
            i < cutoff + 1 && continue

            push!(energy_site, conf.energy / conf.dimension^2)
            push!(E2_site, (conf.energy / conf.dimension^2)^2)
            push!(polar_site, conf.polar / conf.dimension^2)
            push!(abs_polar, abs(conf.polar))
            push!(M2, conf.polar^2)
            push!(M4, conf.polar^4)
        end

        new(energy_site, E2_site, polar_site, abs_polar, M2, M4)
    end
end

struct avgs
    energy_site::Float64
    E2_site::Float64
    polar_site::Float64
    abs_polar::Float64
    M2::Float64
    M4::Float64

    function avgs(macrostate::macrostate)
        energy_site = mean(macrostate.energy_site)
        E2_site = mean(macrostate.E2_site)
        polar_site = mean(macrostate.polar_site)
        abs_polar = mean(macrostate.abs_polar)
        M2 = mean(macrostate.M2)
        M4 = mean(macrostate.M4)
        new(energy_site, E2_site, polar_site, abs_polar, M2, M4)
    end

end
