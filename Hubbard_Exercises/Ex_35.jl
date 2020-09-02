using Plots
import LinearAlgebra.dot

ϵ1 = 2
ϵ0 = 1
t = 1
ax = [1, 0]
ay = [0, 1]

E0(k) = ϵ0
E1(k) =
    (ϵ1 + ϵ0) / 2 +
    1 / 2 * sqrt(
        (ϵ1 + ϵ0)^2 -
        4(ϵ1 * ϵ0 - t^2 * (4 + 2cos(dot(k, ax)) + 2cos(dot(k, ay))))
    )
E2(k) =
    (ϵ1 + ϵ0) / 2 -
    1 / 2 * sqrt(
        (ϵ1 + ϵ0)^2 -
        4(ϵ1 * ϵ0 - t^2 * (4 + 2cos(dot(k, ax)) + 2cos(dot(k, ay))))
    )


points = 1000;
kpoints = [];
num_pts = points^2;

for kx in range(-2π, 2π, length = points)
    for ky in range(-2π, 2π, length = points)
        push!(kpoints, [kx, ky])
    end
end

Energy = []
x = 1
Emax = -Inf
Emin = Inf
for k in kpoints
    global Emax, Emin
    #e0 = E0(k)
    e1 = E1(k)
    e2 = E2(k)
    #Emax < e0 && (Emax = e0)
    Emax < e1 && (Emax = e1)
    Emax < e2 && (Emax = e2)
    #e0 < Emin && (Emin = e0)
    e1 < Emin && (Emin = e1)
    e2 < Emin && (Emin = e2)
    #push!(Energy, e0)
    push!(Energy, e1)
    push!(Energy, e2)
end

println("max is $Emax")
println("min is $Emin")
bin_num = 1000;
bin_size = (Emax - Emin) / bin_num;

hist = zeros(bin_num)
for E in Energy
    index = Int(floor((E - Emin) / bin_size + 1))
    E == Emax && (index = bin_num)
    hist[index] += 1 / (2 * points^2)
end

plot(
    range(Emin, Emax, length = bin_num),
    hist,
    xlabel = "E",
    ylabel = "Density",
    label = "N(E)",
)
