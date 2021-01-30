include("Ising_Driver_Wolff.jl")
using Statistics
using Plots

println("get started")

Temp = (2.1, 2.2, 2.3, 2.4)
sweep = 10^5
Result = driver(12)
println(Result[2])

function error_bunching(Result, bin_size)
    errors = []
    for n = 1:4
        energy_sites = Result[8][n].energy_site
        bin_num = Int(floor(length(energy_sites) / bin_size))
        energy_site_bins = []
        for i = 1:bin_num
            push!(energy_site_bins, mean(energy_sites[((i-1)*bin_size+1):i*bin_size]))
        end
        push!(errors, std(energy_site_bins) / sqrt(bin_num))
    end
    return errors
end


function bunching_plot()
    errors = []
    bin_sizes = [
        2
        4
        8
        16
        32
        64
        128
        256
        512
        1024
        2048
        4096
        8192
        16384
        32768
        65536
    ]
    for bin_size in bin_sizes
        push!(errors, error_bunching(Result, bin_size))
    end
    plot(dpi = 800, title = "error with different bin sizes at different Temperatures")
    for n = 1:4
        errorTn = []
        for m = 1:length(errors)
            push!(errorTn, errors[m][n])
        end
        plot!(
            bin_sizes,
            errorTn,
            label = "T = $(Temp[n])",
            xlabel = "bin size",
            ylabel = "error",
        )

        println(errors)
    end
    savefig("bunching_error_allT.svg")

end

bunching_plot()





#=
struct result
    E::Array
    C::Array
    Binder::Array
    Chi::Array
    function result()
        new([], [], [], [])
    end
end

function Error()
    Temp = range(2.1, 2.4, length = 4)
    println("get started")
    for i = 1:4
        Result = result()
        for j = 1:20
            run = driver(12)
            push!(Result.E, run[2][i])
            push!(Result.C, run[6][i])
            push!(Result.Binder, run[4][i])
            push!(Result.Chi, run[5][i])
        end
        println("T = $(Temp[i])")
        println("T = $(Temp[i]), E = $(mean(Result.E)) +/- $(std(Result.E))")
        println("T = $(Temp[i]), C = $(mean(Result.C)) +/-  $(std(Result.C))")
        println("T = $(Temp[i]), Binder = $(mean(Result.Binder)) +/- $(std(Result.Binder))")
        println("T = $(Temp[i]), Chi = $(mean(Result.Chi)) +/- $(std(Result.Chi))")
    end


end
Error()=#
