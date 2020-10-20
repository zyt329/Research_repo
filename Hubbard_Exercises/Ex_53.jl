using Plots
using LinearAlgebra

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

function Hamiltonian(u, t = 1, T = 1, μ = 0; N::Int64 = 2)
    eye = Matrix{Float64}(I, 4, 4)

    C_up_site = [0 0 0 0; 1 0 0 0; 0 0 0 0; 0 0 1 0]
    D_up_site = [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0]
    C_down_site = [0 0 0 0; 0 0 0 0; 1 0 0 0; 0 1 0 0]
    D_down_site = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
    #=
        C_up_sys = Array{Array{Union{Float64,Int64},2},2}(undef, N, N)
        D_up_sys = Array{Array{Union{Float64,Int64},2},2}(undef, N, N)
        C_down_sys = Array{Array{Union{Float64,Int64},2},2}(undef, N, N)
        D_down_sys = Array{Array{Union{Float64,Int64},2},2}(undef, N, N)
    =#
    C_up_sys = Operater_sys(N, C_up_site)
    D_up_sys = Operater_sys(N, D_up_site)
    C_down_sys = Operater_sys(N, C_down_site)
    D_down_sys = Operater_sys(N, D_down_site)
    @debug println(size(C_up_sys[1, 1]))
    H = zeros(4^(N^2), 4^(N^2))
    for i = 1:N
        for j = 1:N
            H +=
                -t * 1 / 2 * (
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
                )

            H +=
                u * (
                    (
                        C_up_sys[i, j] * D_up_sys[i, j] -
                        0.5 * Matrix{Float64}(I, 4^(N^2), 4^(N^2))
                    ) * (
                        C_down_sys[i, j] * D_down_sys[i, j] -
                        0.5 * Matrix{Float64}(I, 4^(N^2), 4^(N^2))
                    )
                )
            H +=
                -μ * (C_up_sys[i, j] * D_up_sys[i, j] + C_down_sys[i, j] * D_down_sys[i, j])
        end
    end
    return H
end

function M2_gen(N::Int64)
    eye = Matrix{Float64}(I, 4, 4)
    M2_site = [0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0]
    blocks = [M2_site]
    for i = 1:(N^2-1)
        push!(blocks, eye)
    end
    return kron(blocks...)
end

function avg(u, t, T, μ; N::Int64)
    M2 = M2_gen(N)
    @assert(
        Hermitian(Hamiltonian(u, t, T, μ, N = N)) == Hamiltonian(u, t, T, μ, N = N),
        "Hamiltonian is not Hermitian."
    )
    H = Hermitian(Hamiltonian(u, t, T, μ, N = N))
    eigH = eigen(H)
    #println(eigH.values)
    ρ = diagm(0 => exp.(-eigH.values / T))
    E = diagm(0 => eigH.values)
    ρ = ρ / tr(ρ)
    @assert(Diagonal(ρ) == ρ, "ρ is not diagonal")
    #=@debug
    println(tr(ρ))
    for element in diag(ρ)
        element > 1 && println("\n", element)
    end=#
    #println("\n", diag(ρ .* E))
    #println("\n", diag(E))
    E_avg = tr(ρ .* E)
    M2_avg = tr(ρ * eigH.vectors' * M2 * eigH.vectors)
    return (E_avg, M2_avg)
end

println(avg(100, 1, 1, 0, N = 2)[1])

p1 = plot(
    1:256,
    eigen(Hermitian(Hamiltonian(1, 1, 100, 0, N = 2))).values,
    ylabel = "Energy",
    label = "Energy Levels(U=1,t=1,T=1,μ=0)",
    title = "Energy Levels",
)
plot!(
    1:256,
    eigen(Hermitian(Hamiltonian(5, 1, 100, 0, N = 2))).values,
    label = "Energy Levels(U=5,t=1,T=1,μ=0)",
)
plot!(
    1:256,
    eigen(Hermitian(Hamiltonian(10, 1, 100, 0, N = 2))).values,
    label = "Energy Levels(U=10,t=1,T=1,μ=0)",
)

M2u(u) = avg(u, 1, 1, 0, N = 2)[2]
p2 = plot(M2u, 0, 100, xlabel = "U", ylabel = "M^2", label = "M^2(t=1,T=1)", title = "M^2")

E(T) = avg(1, 1, T, 0, N = 2)[1]
p3 = plot(E, 0.1, 100, xlabel = "T", ylabel = "E", label = "E(t=1,U=1)", title = "E")

plot(p1, p2, p3, layout = (3, 1))
