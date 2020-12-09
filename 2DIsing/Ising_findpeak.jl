using JLD
using Dates

Result = load("Result_Thu_03_Dec_2020_01_35_45.jld")["Result"]
#push!(Result,load("Result4_Wed_02_Dec_2020_16_03_35.jld")["Result"])

"""
Finding the maximum(peak) of Array A by averaging over its "avgnum" of its neighbours on both of its sides. (Ignoring the boundaries)

Output:
(maxval, index)
maxval: Maximum value found at Array A.
index: The index on which we found the maximum of A.
"""
function findpeak(A::Array, neighbour::Int = 0)
    @assert(
        neighbour < length(A) / 2,
        "You are averaging over too many neighbours! The neighbours that are averaged over should be less than half the length of the Array."
    )
    Δi = neighbour
    B = []
    for n = (1+Δi):(length(A)-Δi)
        avg = sum(A[n-Δi:n+Δi]) / (2Δi + 1)
        push!(B, avg)
    end
    maxim = findmax(B)
    return (A[maxim[2]+Δi], maxim[2] + Δi)
end

"""
Finding the maximum of χ and C from the result got from running the driver.

Input:
Result:Result.jld

Output:
peaks: It's an array of two entries. Each entry of the array is another array that represents the peak information. [[...],[...]]
The first array represents peaks of χ in different sizes of lattices. The second array represents peaks of C in different sizes of lattices.
Each entry of the original array is an array of tuples. Each tuple has two arguments. The first argument gives the maximum value of χ/C, the second argument gives the temperature of the peak.
For example, it is like this:[[(peak_val_of_χ_8*8,Temp_for_peak_of_χ_8*8),(...12*12),(...16*16)],[(...of_C_8*8),(...12*12),(...16*16)]]
"""
function findmaxResult(Result)
    println("   ")
    peaks = [[], []]
    for i = 1:4
        maxim = findpeak(Result[i][5], 5)
        push!(peaks[1], (maxim[1], Result[i][1][maxim[2]]))
        println("maximum of χ_site_$(Result[i][7]) is $(maxim[1]), with a temperature of $(Result[i][1][maxim[2]])")
        maxim = findpeak(Result[i][6], 5)
        push!(peaks[2], (maxim[1], Result[i][1][maxim[2]]))
        println("maximum of C_site_$(Result[i][7]) is $(maxim[1]), with a temperature of $(Result[i][1][maxim[2]])")
    end
    return peaks
end

peaks = findmaxResult(Result)
my_time = Dates.now()
save("peaks_$(Dates.format(my_time, "e_dd_u_yyyy_HH_MM_SS")).jld", "peaks", peaks)

#findmaxResult(Result)
(
    (2.543859649122807, 1.7634769316172738),
    (2.4210526315789473, 2.9311296858223046),
    (2.433333333333333, 3.4716471650590264),
    (2.3473684210526318, 3.1245134180556198),
    (2.3964912280701753, 5.56867382635552),
    (2.3350877192982455, 3.6540187034766736),
)

[
    [
        (1.7634769316172738, 2.543859649122807),
        (3.4716471650590264, 2.433333333333333),
        (5.56867382635552, 2.3964912280701753),
        (13.67053374738642, 2.3719298245614033),
    ],
    [
        (2.9311296858223046, 2.4210526315789473),
        (3.1245134180556198, 2.3473684210526318),
        (3.6540187034766736, 2.3350877192982455),
        (3.8776559405177538, 2.3105263157894735),
    ],
]
