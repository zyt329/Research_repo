using Plots

t = 1
v = 0.1

E1(k) = -t * cos(k) + sqrt(t^2 * cos(k)^2 + v^2)
E2(k) = -t * cos(k) - sqrt(t^2 * cos(k)^2 + v^2)

p1 = plot(E1, -2π, 2π)
plot!(E2, -2π, 2π)

points = 1000000;
t = 1;
kpoints = [];
num_pts = points^2;

for k in range(-2π, 2π, length = points)
        push!(kpoints, k)
end

Energy = []
x = 1
Emax = -Inf
Emin = Inf
for k in kpoints
    global Emax,Emin
    E0 = E1(k)
    E00 = E2(k)
    Emax < E0 && (Emax = E0)
    Emax < E00 && (Emax = E00)
    E0 < Emin && (Emin = E0)
    E00 < Emin && (Emin = E00)
    push!(Energy, E0)
    push!(Energy, E00)
end

println("max is $Emax")
println("min is $Emin")
bin_num = 1000;
bin_size = (Emax - Emin) / bin_num ;

hist = zeros(bin_num)
for E in Energy
    index = Int( floor( (E - Emin) / bin_size + 1 ) )
    E == Emax && (index = bin_num)
    hist[index] += 1 / (2 * points)
end

p2 = plot(
    range(Emin, Emax, length = bin_num),
    hist,
    xlabel = "E",
    ylabel = "Density",
    label = "N(E)",
)

plot(p1,p2,layout=(2,1))
