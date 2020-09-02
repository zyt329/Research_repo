using Plots

points = 10000;
t = 1;
kpoints = [];
num_pts = points^2;
bin_num = 10000;
bin_size = 8t / bin_num;

for kx in range(-2π, 2π, length = points)
    for ky in range(-2π, 2π, length = points)
        push!(kpoints, (kx, ky))
    end
end

E(kpoints) = -2t * (cos(kpoints[1]) + cos(kpoints[2]))

hist = zeros(bin_num)
for k in kpoints
    index = Int(floor( (E(k) + 4t)/bin_size ) + 1)
    hist[index] += 1 / points^2
end

plot(
    range(-4t, 4t, length = bin_num),
    hist,
    xlabel = "E",
    ylabel = "Density",
    label = "N(E)",
)
