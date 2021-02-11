using Plots
import LinearAlgebra.dot

points = 2000;
t = 1;
kpoints = [];
num_pts = points^2;

for kx in range(-2π, 2π, length = points)
    for ky in range(-2π, 2π, length = points)
        push!(kpoints, [kx, ky])
    end
end

#ax = [1, 0]; ay = [1 / 2, sqrt(3) / 2]
E(k) = 2 * cos(k[1] - k[2]) + 2 * cos(k[1]) + 2 * cos(k[2])

Energy = []
x = 1
Emax = -Inf
Emin = Inf
for k in kpoints
    global Emax,Emin
    E0 = E(k)
    Emax < E0 && (Emax = E0)
    E0 < Emin && (Emin = E0)
    push!(Energy, E0)
end

println("max is $Emax")
println("min is $Emin")
bin_num = 1000;
bin_size = (Emax - Emin) / bin_num ;

hist = zeros(bin_num)
for E in Energy
    index = Int( floor( (E - Emin) / bin_size + 1 ) )
    E == Emax && (index = bin_num)
    hist[index] += 1 / points^2
end

plot(
    range(Emin, Emax, length = bin_num),
    hist,
    xlabel = "E",
    ylabel = "Density",
    label = "N(E)",
)

println("Total energy is ",sum(range(Emin, Emax, length = bin_num).*hist))
