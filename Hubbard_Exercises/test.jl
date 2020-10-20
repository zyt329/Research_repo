using Plots
using LinearAlgebra

t = 1

eye = Matrix{Float64}(I, 4, 4)

C_up_site = [0 0 0 0; 1 0 0 0; 0 0 0 0; 0 0 1 0]
D_up_site = [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0]
C_down_site = [0 0 0 0; 0 0 0 0; 1 0 0 0; 0 1 0 0]
D_down_site = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]

function Operater_sys(N::Int, C_site::Array)
    Op_sys = Array{Array,2}(undef, N, N)
    for i = 1:N
        for j = 1:N
            blocks = []
            for l = 1:N^2
                if l == ((i - 1) * N + j)
                    push!(blocks, C_site)
                    continue
                end
                push!(blocks, eye)
            end
            Op_sys[i, j] = kron(blocks...)
        end
    end
    return Op_sys
end

#=
println(sum(C_up_sys[2,2]))
println(typeof(C_up_sys[1,1]))
println(size(C_up_sys[1,1]*D_up_sys[1,1]))
=#
function Hamiltonian(N)
    C_up_sys = Operater_sys(N, C_up_site)
    D_up_sys = Operater_sys(N, D_up_site)
    C_down_sys = Operater_sys(N, C_down_site)
    D_down_sys = Operater_sys(N, D_down_site)

    H = zeros(4^(N^2), 4^(N^2))
    for i = 1:N
        for j = 1:N
        #=    H += 1/2 * (
                C_up_sys[i, j] * (
                    D_up_sys[mod1(i - 1, N), j] +
                    D_up_sys[mod1(i + 1, N), j] +
                    D_up_sys[i, mod1(j - 1, N)] +
                    D_up_sys[i, mod1(j + 1, N)]
                ) +
                C_down_sys[i, j] * (
                    D_down_sys[mod1(i - 1, N), j] +
                    D_down_sys[mod1(i + 1, N), j] +
                    D_down_sys[i, mod1(j - 1, N)] +
                    D_down_sys[i, mod1(j + 1, N)]
                )
            )=#

            H +=
                 0.5 * (
                    (
                        C_up_sys[i, j] * D_up_sys[i, j] -
                        0.5 * Matrix{Float64}(I, 4^(N^2), 4^(N^2))
                    ) * (
                        C_down_sys[i, j] * D_down_sys[i, j] -
                        0.5 * Matrix{Float64}(I, 4^(N^2), 4^(N^2))
                    )
                )
            H +=
                0.1 * (C_up_sys[i, j] * D_up_sys[i, j] + C_down_sys[i, j] * D_down_sys[i, j])
        end
    end
    return H
end

println(typeof(Hamiltonian(2)))


@assert(Diagonal(Hamiltonian(2)) == Hamiltonian(2), "Hamiltonian is not diagonal.")

println("Check done.")
